#include "spin2014/SlqpIterationStore.h"
#include "manufacturing/SlicerValidation.h"
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <climits>
#include <fcntl.h>
#include <filesystem>
#include <limits>
#include <sys/stat.h>
#include <type_traits>
#include <unistd.h>

namespace spin2014 {
namespace {
constexpr std::uint64_t recordLimit=64ULL*1024*1024;
[[noreturn]] void error(const std::string& message){throw SlqpStorageError("SLQP iteration storage: "+message);}
void appendU64(std::string& b,std::uint64_t v){for(int i=0;i<8;++i)b.push_back(char(v>>(8*i)));}
std::uint64_t loadU64(const char* p){std::uint64_t n=0;for(int i=0;i<8;++i)n|=std::uint64_t(static_cast<unsigned char>(p[i]))<<(8*i);return n;}
void readExactAt(int descriptor,char* data,std::size_t n,std::uint64_t offset) {
    std::size_t done=0;
    while(done<n) {
        const auto got=::pread(descriptor,data+done,n-done,static_cast<off_t>(offset+done));
        if(got<0 && errno==EINTR)continue;
        if(got<=0)error("truncated/read failure");
        done+=static_cast<std::size_t>(got);
    }
}
std::uint64_t readU64At(int descriptor,std::uint64_t offset){char bytes[8];readExactAt(descriptor,bytes,8,offset);return loadU64(bytes);}
bool sameStatSnapshot(const struct stat& a,const struct stat& b) {
    return S_ISREG(a.st_mode) && S_ISREG(b.st_mode) && a.st_dev==b.st_dev && a.st_ino==b.st_ino &&
        a.st_size==b.st_size && a.st_mtim.tv_sec==b.st_mtim.tv_sec && a.st_mtim.tv_nsec==b.st_mtim.tv_nsec &&
        a.st_ctim.tv_sec==b.st_ctim.tv_sec && a.st_ctim.tv_nsec==b.st_ctim.tv_nsec;
}
bool samePathIdentity(const struct stat& a,const struct stat& b) {
    return S_ISREG(b.st_mode) && a.st_dev==b.st_dev && a.st_ino==b.st_ino;
}

// One field traversal for both encoder and decoder prevents asymmetric fields.
template<class A,class T> void fields(A& a,T& t){
    a(t.iteration);a(t.objective);a(t.equalityInfinityNorm);a(t.projectedGradientInfinityNorm);
    a(t.kktInfinityNorm);a(t.meritValue);a(t.acceptedStep);a(t.lineSearchEvaluations);
    a(t.activeVariables);a(t.lbfgsHistorySize);a(t.usedFirstOrderFallback);a(t.phi);
    a(t.kkt);a(t.state);a(t.gradient);a(t.constraints);a(t.constraintJacobian);
    a(t.activeMembership);a(t.membershipHash);a(t.historyResetReason);a(t.historyBefore);a(t.historyAfter);
    a(t.curvaturePairAccepted);a(t.curvaturePairRejected);a(t.linearizedConstraintResidual);
    a(t.restoration);a(t.objectiveTrials);a(t.finalState);a(t.metricTangent);
    a(t.metricTangentSelected);a(t.metricTangentAccepted);a(t.metricTangentTrialCount);a(t.metricTangentFallbackAttempted);
}
template<class A,class T> void kktFields(A& a,T& t){a(t.valid);a(t.status);a(t.lambda);a(t.alpha);a(t.eta);a(t.equalityPrimal);a(t.boundPrimal);a(t.stationarity);a(t.dualSign);a(t.complementarity);a(t.full);}
template<class A,class T> void trialFields(A& a,T& t){a(t.step);a(t.constraintL1);a(t.valid);a(t.accepted);a(t.reason);}
template<class A,class T> void restorationFields(A& a,T& t){a(t.attempted);a(t.accepted);a(t.status);a(t.reason);a(t.direction);a(t.predictedResidual);a(t.trials);}
template<class A,class T> void tangentFields(A& a,T& t){a(t.attempted);a(t.checked);a(t.status);a(t.direction);a(t.workingSet);a(t.multipliers);a(t.equalityResidual);a(t.stationarity);a(t.trace);}
template<class A,class T> void stepFields(A& a,T& t){a(t.workingSet);a(t.action);a(t.changedIndex);a(t.equalityResidual);a(t.stationarity);}
#define STRUCT_FIELDS(QUAL) \
 void operator()(QUAL KktResidual& t){kktFields(*this,t);} \
 void operator()(QUAL RestorationTrial& t){trialFields(*this,t);} \
 void operator()(QUAL RestorationEvidence& t){restorationFields(*this,t);} \
 void operator()(QUAL TangentQpResult& t){tangentFields(*this,t);} \
 void operator()(QUAL TangentQpStep& t){stepFields(*this,t);}
struct Encoder {
    std::string bytes;
    void room(std::uint64_t n){if(n>recordLimit-bytes.size())error("record exceeds 64 MiB format bound");}
    void operator()(std::uint64_t n){room(8);appendU64(bytes,n);}
    void operator()(int n){(*this)(std::uint64_t(std::int64_t(n)));}
    void operator()(bool b){room(1);bytes.push_back(b?1:0);}
    void operator()(unsigned char b){room(1);bytes.push_back(char(b));}
    void operator()(double x){std::uint64_t bits;static_assert(sizeof(x)==sizeof(bits));std::memcpy(&bits,&x,8);(*this)(bits);}
    void operator()(const std::string& s){(*this)(std::uint64_t(s.size()));room(s.size());bytes.append(s);}
    template<class T> void operator()(const std::vector<T>& v){(*this)(std::uint64_t(v.size()));for(const auto& t:v)(*this)(t);}
    template<class T> void matrix(const T& v){(*this)(std::uint64_t(v.rows()));(*this)(std::uint64_t(v.cols()));for(Eigen::Index i=0;i<v.size();++i)(*this)(v.data()[i]);}
    void operator()(const Eigen::VectorXd& v){matrix(v);}
    void operator()(const Eigen::MatrixXd& v){matrix(v);}
    STRUCT_FIELDS(const)
};
struct Decoder {
    const std::string& bytes;std::size_t pos=0;std::uint64_t allocations=0;
    void need(std::uint64_t n){if(n>bytes.size()-pos)error("invalid record length");}
    void allocate(std::uint64_t n){if(n>2*recordLimit-allocations)error("record allocation bound exceeded");allocations+=n;}
    void operator()(std::uint64_t& n){need(8);n=loadU64(bytes.data()+pos);pos+=8;}
    void operator()(int& n){std::uint64_t v;(*this)(v);std::int64_t signedValue;std::memcpy(&signedValue,&v,8);if(signedValue<INT_MIN||signedValue>INT_MAX)error("invalid integer");n=int(signedValue);}
    void operator()(bool& b){need(1);auto c=bytes[pos++];if(c!=0&&c!=1)error("invalid boolean");b=c;}
    void operator()(unsigned char& b){need(1);b=static_cast<unsigned char>(bytes[pos++]);}
    void operator()(double& x){std::uint64_t bits;(*this)(bits);std::memcpy(&x,&bits,8);}
    void operator()(std::string& s){std::uint64_t n;(*this)(n);need(n);allocate(n);s.assign(bytes.data()+pos,n);pos+=n;}
    template<class T> void operator()(std::vector<T>& v){std::uint64_t n;(*this)(n);need(n);if(n>2*recordLimit/sizeof(T))error("invalid array size");allocate(n*sizeof(T));v.resize(n);for(auto& t:v)(*this)(t);}
    template<class T> void matrix(T& v,bool vector){std::uint64_t rows,cols;(*this)(rows);(*this)(cols);
        if(rows>recordLimit/8||cols>recordLimit/8||(vector&&cols!=1)|| (cols&&rows>recordLimit/8/cols))error("invalid matrix dimensions");
        need(rows*cols*8);allocate(rows*cols*8);v.resize(rows,cols);for(Eigen::Index i=0;i<v.size();++i)(*this)(v.data()[i]);}
    void operator()(Eigen::VectorXd& v){matrix(v,true);}
    void operator()(Eigen::MatrixXd& v){matrix(v,false);}
    STRUCT_FIELDS()
};
#undef STRUCT_FIELDS
}
std::string encodeSlqpIteration(const SlqpIteration& t){Encoder e;fields(e,t);return std::move(e.bytes);}
SlqpIteration decodeSlqpIteration(const std::string& bytes) {
    if(bytes.size()>recordLimit)error("record exceeds format bound");
    SlqpIteration t;Decoder d{bytes};fields(d,t);
    if(d.pos!=bytes.size())error("trailing record bytes");return t;
}
std::shared_ptr<SlqpIterationStore> SlqpIterationStore::openSealed(
    const std::string& path,std::uint64_t maximumBytes) {
    return openSealedAndVisit(path,maximumBytes,{});
}
std::shared_ptr<SlqpIterationStore> SlqpIterationStore::openSealedAndVisit(
    const std::string& path,std::uint64_t maximumBytes,
    const std::function<void(const SlqpIteration&)>& visit) {
    const auto anchored=std::filesystem::absolute(path);
    const int descriptor=::open(anchored.c_str(),O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if(descriptor<0)error("cannot open read-only history without following symlinks: "+anchored.string());
    int ownedDescriptor=descriptor;
    try {
        struct stat before{},pathBefore{};
        if(::fstat(descriptor,&before)!=0 || !S_ISREG(before.st_mode) || before.st_size<0)
            error("read-only history requires a regular file");
        if(static_cast<std::uint64_t>(before.st_size)>maximumBytes)
            error("read-only history exceeds remaining byte budget");
        if(::lstat(anchored.c_str(),&pathBefore)!=0 || !samePathIdentity(before,pathBefore))
            error("history pathname changed before reopening");
        SlqpIterationReader reader(descriptor);SlqpIteration t;std::size_t count=0;
        if(reader.device_!=static_cast<std::uint64_t>(before.st_dev) || reader.inode_!=static_cast<std::uint64_t>(before.st_ino) ||
           reader.fileSize_!=static_cast<std::uint64_t>(before.st_size))
            error("history changed while its reader was initialized");
        while(reader.next(t)) {++count;if(visit)visit(t);}
        struct stat after{},pathAfter{};
        if(::fstat(descriptor,&after)!=0 || !sameStatSnapshot(before,after) || !reader.snapshotUnchanged())
            error("history changed while being reopened");
        if(::lstat(anchored.c_str(),&pathAfter)!=0 || !samePathIdentity(before,pathAfter) ||
           pathAfter.st_size!=before.st_size)
            error("history pathname changed while being reopened");
        auto* raw=new SlqpIterationStore(anchored.string(),maximumBytes,ReadOnlyTag{},descriptor,reader);
        ownedDescriptor=-1;
        auto out=std::shared_ptr<SlqpIterationStore>(raw);
        out->count_=count;out->bytes_=static_cast<std::uint64_t>(before.st_size);out->sealed_=true;
        return out;
    } catch(...) {
        if(ownedDescriptor>=0)::close(ownedDescriptor);
        throw;
    }
}
std::unique_ptr<SlqpIterationReader> SlqpIterationStore::openReader() const {
    if(!sealed_ || readOnlyDescriptor_<0)error("read-only history has no validated file descriptor");
    const int copy=::fcntl(readOnlyDescriptor_,F_DUPFD_CLOEXEC,0);
    if(copy<0)error("cannot duplicate validated history descriptor");
    struct stat current{};
    if(::fstat(copy,&current)!=0 || !S_ISREG(current.st_mode) ||
       static_cast<std::uint64_t>(current.st_dev)!=device_ || static_cast<std::uint64_t>(current.st_ino)!=inode_ ||
       static_cast<std::uint64_t>(current.st_size)!=bytes_ || current.st_mtim.tv_sec!=mtimeSeconds_ ||
       current.st_mtim.tv_nsec!=mtimeNanoseconds_ || current.st_ctim.tv_sec!=ctimeSeconds_ ||
       current.st_ctim.tv_nsec!=ctimeNanoseconds_) {
        ::close(copy);error("validated history changed in place");
    }
    struct stat named{};
    if(::lstat(path_.c_str(),&named)!=0 || !S_ISREG(named.st_mode) || named.st_dev!=current.st_dev ||
       named.st_ino!=current.st_ino || named.st_size!=current.st_size) {
        ::close(copy);error("validated history pathname was replaced");
    }
    std::unique_ptr<SlqpIterationReader> reader;
    try {reader.reset(new SlqpIterationReader(copy,SlqpIterationReader::AdoptDescriptorTag{}));}
    catch(...) {::close(copy);throw;}
    if(reader->device_!=device_ || reader->inode_!=inode_ || reader->fileSize_!=bytes_ ||
       reader->mtimeSeconds_!=mtimeSeconds_ || reader->mtimeNanoseconds_!=mtimeNanoseconds_ ||
       reader->ctimeSeconds_!=ctimeSeconds_ || reader->ctimeNanoseconds_!=ctimeNanoseconds_)
        error("validated history changed while opening a reader");
    return reader;
}
void SlqpIterationStore::validateUnchangedSnapshot() const {
    // initialize() checks the pinned descriptor, pathname identity, header,
    // footer and captured stat snapshot; intentionally do not call next().
    // The caller must already hold a digest from a full validated traversal.
    (void)openReader();
}
SlqpIterationStore::SlqpIterationStore(std::string path,std::uint64_t maximumBytes):path_(std::move(path)),maximum_(maximumBytes){
    std::error_code pathError;
    const auto anchored=std::filesystem::absolute(path_,pathError);
    if(pathError)error("cannot anchor artifact path: "+pathError.message());
    path_=anchored.string();
    const int fd=::open(path_.c_str(),O_RDWR|O_CREAT|O_EXCL|O_CLOEXEC,0600);
    if(fd<0)error("exclusive creation failed: "+path_);
    output_=::fdopen(fd,"wb");if(!output_){::close(fd);error("opening writer failed: "+path_);}
    // Disable stdio buffering so write failures and bytesWritten are observable.
    std::setvbuf(output_,nullptr,_IONBF,0);
    try{std::string header="SLQPTRC1";appendU64(header,1);appendU64(header,recordLimit);write(header.data(),header.size());}
    catch(...){std::fclose(output_);output_=nullptr;throw;}
}
SlqpIterationStore::SlqpIterationStore(std::string path,std::uint64_t maximumBytes,ReadOnlyTag,
    int descriptor,const SlqpIterationReader& snapshot)
    :path_(std::move(path)),readOnlyDescriptor_(descriptor),maximum_(maximumBytes),
     device_(snapshot.device_),inode_(snapshot.inode_),mtimeSeconds_(snapshot.mtimeSeconds_),
     mtimeNanoseconds_(snapshot.mtimeNanoseconds_),ctimeSeconds_(snapshot.ctimeSeconds_),
     ctimeNanoseconds_(snapshot.ctimeNanoseconds_) {}
SlqpIterationStore::~SlqpIterationStore(){if(output_)std::fclose(output_);if(readOnlyDescriptor_>=0)::close(readOnlyDescriptor_);}
void SlqpIterationStore::write(const void* p,std::size_t n){
    if(failed_||sealed_)error("writer failed or sealed: "+path_);
    if(n>maximum_-bytes_){failed_=true;error("byte resource limit: "+path_);}
    const auto written=std::fwrite(p,1,n,output_);bytes_+=written;
    if(written!=n){failed_=true;error("write failure; partial evidence preserved: "+path_);}
}
void SlqpIterationStore::flushPending(){if(!pending_)return;
    try{auto payload=encodeSlqpIteration(*pending_);auto digest=manufacturing::computeContentSha256(payload);
        std::string size;appendU64(size,payload.size());write(size.data(),size.size());write(digest.data(),digest.size());write(payload.data(),payload.size());pending_.reset();
    }catch(...){failed_=true;throw;}
}
void SlqpIterationStore::append(SlqpIteration t){if(sealed_||failed_)error("append to sealed/failed store");flushPending();pending_=std::move(t);++count_;peak_=1;}
SlqpIteration& SlqpIterationStore::back(){if(!pending_||sealed_||failed_)error("no mutable pending record");return *pending_;}
void SlqpIterationStore::seal(){if(sealed_||failed_)error("seal on sealed/failed store");flushPending();
    if(std::fflush(output_)!=0||::fsync(::fileno(output_))!=0){failed_=true;error("data sync failure: "+path_);}
    const auto footerOffset=bytes_;
    std::string footer="SLQPEND1";appendU64(footer,count_);appendU64(footer,bytes_);write(footer.data(),footer.size());
    if(std::fflush(output_)!=0||::fsync(::fileno(output_))!=0){
        // A reported sync failure must not leave a readable completion marker.
        // Preserve all data records, removing only this failed seal attempt.
        const bool invalidated=::ftruncate(::fileno(output_),footerOffset)==0;
        failed_=true;error(std::string("seal sync failure; ")+(invalidated?"seal invalidated: ":"seal invalidation also failed: ")+path_);
    }
    const int pinned=::fcntl(::fileno(output_),F_DUPFD_CLOEXEC,0);
    if(pinned<0) {
        const bool invalidated=::ftruncate(::fileno(output_),footerOffset)==0;
        failed_=true;error(std::string("descriptor pin failed; ")+(invalidated?"seal invalidated: ":"seal invalidation also failed: ")+path_);
    }
    struct stat metadata{};
    if(::fstat(pinned,&metadata)!=0 || !S_ISREG(metadata.st_mode) || metadata.st_size<0) {
        ::close(pinned);const bool invalidated=::ftruncate(::fileno(output_),footerOffset)==0;
        failed_=true;error(std::string("sealed descriptor validation failed; ")+(invalidated?"seal invalidated: ":"seal invalidation also failed: ")+path_);
    }
    const auto closed=std::fclose(output_);output_=nullptr;
    if(closed!=0){::close(pinned);failed_=true;error("close failure: "+path_);}
    readOnlyDescriptor_=pinned;device_=metadata.st_dev;inode_=metadata.st_ino;
    mtimeSeconds_=metadata.st_mtim.tv_sec;mtimeNanoseconds_=metadata.st_mtim.tv_nsec;
    ctimeSeconds_=metadata.st_ctim.tv_sec;ctimeNanoseconds_=metadata.st_ctim.tv_nsec;sealed_=true;
}
SlqpIterationReader::SlqpIterationReader(const std::string& path) {
    descriptor_=::open(path.c_str(),O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if(descriptor_<0)error("opening reader failed: "+path);
    try{initialize();}catch(...){::close(descriptor_);descriptor_=-1;throw;}
}
SlqpIterationReader::SlqpIterationReader(int descriptor)
    :descriptor_(::fcntl(descriptor,F_DUPFD_CLOEXEC,0)) {
    if(descriptor_<0)error("duplicating reader descriptor failed");
    try{initialize();}catch(...){::close(descriptor_);descriptor_=-1;throw;}
}
SlqpIterationReader::SlqpIterationReader(int descriptor,AdoptDescriptorTag):descriptor_(descriptor) {
    if(descriptor_<0)error("duplicating reader descriptor failed");
    initialize();
}
SlqpIterationReader::~SlqpIterationReader(){if(descriptor_>=0)::close(descriptor_);}
void SlqpIterationReader::initialize() {
    struct stat metadata{};if(::fstat(descriptor_,&metadata)!=0 || !S_ISREG(metadata.st_mode) || metadata.st_size<0)
        error("reader requires a regular file");
    fileSize_=static_cast<std::uint64_t>(metadata.st_size);device_=static_cast<std::uint64_t>(metadata.st_dev);inode_=static_cast<std::uint64_t>(metadata.st_ino);
    mtimeSeconds_=metadata.st_mtim.tv_sec;mtimeNanoseconds_=metadata.st_mtim.tv_nsec;
    ctimeSeconds_=metadata.st_ctim.tv_sec;ctimeNanoseconds_=metadata.st_ctim.tv_nsec;
    if(fileSize_<48)error("missing completion");
    char magic[8];readExactAt(descriptor_,magic,8,0);
    if(std::memcmp(magic,"SLQPTRC1",8)||readU64At(descriptor_,8)!=1||readU64At(descriptor_,16)!=recordLimit)
        error("unsupported format/version");
    readExactAt(descriptor_,magic,8,fileSize_-24);if(std::memcmp(magic,"SLQPEND1",8))error("missing completion");
    remaining_=readU64At(descriptor_,fileSize_-16);end_=readU64At(descriptor_,fileSize_-8);
    if(end_!=fileSize_-24||remaining_>(end_-24)/72)error("invalid completion/count");
    cursor_=24;
}
bool SlqpIterationReader::snapshotUnchanged() const {
    struct stat current{};
    return ::fstat(descriptor_,&current)==0 && S_ISREG(current.st_mode) &&
        static_cast<std::uint64_t>(current.st_dev)==device_ && static_cast<std::uint64_t>(current.st_ino)==inode_ &&
        static_cast<std::uint64_t>(current.st_size)==fileSize_ && current.st_mtim.tv_sec==mtimeSeconds_ &&
        current.st_mtim.tv_nsec==mtimeNanoseconds_ && current.st_ctim.tv_sec==ctimeSeconds_ &&
        current.st_ctim.tv_nsec==ctimeNanoseconds_;
}
bool SlqpIterationReader::next(SlqpIteration& t){
    if(cursor_>end_)error("invalid read position");
    if(!remaining_){if(cursor_!=end_)error("completion count mismatch");if(!snapshotUnchanged())error("history changed while reading");return false;}
    if(end_-cursor_<72)error("truncated record frame");
    const auto size=readU64At(descriptor_,cursor_);if(size>recordLimit||size>end_-cursor_-72)error("invalid record size");cursor_+=8;
    std::string digest(64,'\0');readExactAt(descriptor_,digest.data(),64,cursor_);cursor_+=64;
    std::string bytes(size,'\0');readExactAt(descriptor_,bytes.data(),size,cursor_);cursor_+=size;
    if(manufacturing::computeContentSha256(bytes)!=digest)error("record checksum mismatch");
    t=decodeSlqpIteration(bytes);--remaining_;
    if(!remaining_ && (cursor_!=end_ || !snapshotUnchanged()))error("completion mismatch or history changed while reading");
    return true;
}
std::size_t slqpTraceCount(const SlqpLbfgsResult& r){return r.iterationStore?r.iterationStore->count():r.trace.size();}
void forEachSlqpIteration(const SlqpLbfgsResult& r,const std::function<void(const SlqpIteration&)>& visit){
    if(r.iterationStore){auto reader=r.iterationStore->openReader();SlqpIteration t;while(reader->next(t))visit(t);}
    else for(const auto& t:r.trace)visit(t);
}
}
