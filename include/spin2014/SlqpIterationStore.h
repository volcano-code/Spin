#pragma once
#include "spin2014/SlqpLbfgsSolver.h"
#include <cstdio>
#include <fstream>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>

namespace spin2014 {
// This format is a versioned local binary artifact, not the JSONL archive.
// Every field is encoded; per-record lengths and SHA-256 protect readback.
class SlqpStorageError : public std::runtime_error {
public: using std::runtime_error::runtime_error;
};
std::string encodeSlqpIteration(const SlqpIteration&);
SlqpIteration decodeSlqpIteration(const std::string&);
class SlqpIterationReader {
public:
    explicit SlqpIterationReader(const std::string& path);
    // Duplicates an already-open descriptor so the reader is independent of
    // later pathname replacement. The duplicate is owned by this reader.
    explicit SlqpIterationReader(int descriptor);
    ~SlqpIterationReader();
    SlqpIterationReader(const SlqpIterationReader&)=delete;
    SlqpIterationReader& operator=(const SlqpIterationReader&)=delete;
    bool next(SlqpIteration&);
private:
    friend class SlqpIterationStore;
    struct AdoptDescriptorTag {};
    SlqpIterationReader(int descriptor,AdoptDescriptorTag);
    void initialize();
    bool snapshotUnchanged() const;
    int descriptor_=-1;
    std::uint64_t remaining_=0, end_=0, cursor_=24, fileSize_=0;
    std::uint64_t device_=0,inode_=0;
    std::int64_t mtimeSeconds_=0,mtimeNanoseconds_=0,ctimeSeconds_=0,ctimeNanoseconds_=0;
};
class SlqpIterationStore {
public:
    explicit SlqpIterationStore(std::string path, std::uint64_t maximumBytes);
    // No writer is opened. Validate all records/completion within the size
    // cap before returning an immutable (sealed) streaming history handle.
    static std::shared_ptr<SlqpIterationStore> openSealed(
        const std::string& path,std::uint64_t maximumBytes);
    // Validate and visit each record during the same streaming pass. The
    // callback runs only after that record's frame, digest, and payload have
    // been validated; the store is returned only after the complete footer
    // and file-size checks succeed.
    static std::shared_ptr<SlqpIterationStore> openSealedAndVisit(
        const std::string& path,std::uint64_t maximumBytes,
        const std::function<void(const SlqpIteration&)>& visit);
    // Opens another reader on the exact validated inode. It refuses a file
    // modified in place since sealing; later traversal still rechecks every
    // record checksum and file metadata.
    std::unique_ptr<SlqpIterationReader> openReader() const;
    // Rechecks the pinned descriptor and pathname snapshot without traversing
    // record payloads. Safe only when the caller already holds a digest that
    // was computed while validating this exact sealed snapshot.
    void validateUnchangedSnapshot() const;
    ~SlqpIterationStore();
    SlqpIterationStore(const SlqpIterationStore&)=delete;
    SlqpIterationStore& operator=(const SlqpIterationStore&)=delete;
    void append(SlqpIteration);
    SlqpIteration& back();
    void seal();
    std::size_t count() const {return count_;}
    std::size_t retainedRecords() const {return pending_.has_value();}
    std::size_t peakRetainedRecords() const {return peak_;}
    std::uint64_t bytesWritten() const {return bytes_;}
    const std::string& path() const {return path_;}
    bool sealed() const {return sealed_;}
private:
    struct ReadOnlyTag {};
    SlqpIterationStore(std::string path,std::uint64_t maximumBytes,ReadOnlyTag,
        int descriptor,const SlqpIterationReader& snapshot);
    void write(const void*,std::size_t);
    void flushPending();
    std::string path_;
    std::FILE* output_=nullptr;
    int readOnlyDescriptor_=-1;
    std::optional<SlqpIteration> pending_;
    std::uint64_t maximum_,bytes_=0;
    std::size_t count_=0,peak_=0;
    std::uint64_t device_=0,inode_=0;
    std::int64_t mtimeSeconds_=0,mtimeNanoseconds_=0,ctimeSeconds_=0,ctimeNanoseconds_=0;
    bool sealed_=false,failed_=false;
};
std::size_t slqpTraceCount(const SlqpLbfgsResult&);
void forEachSlqpIteration(const SlqpLbfgsResult&,const std::function<void(const SlqpIteration&)>&);
}
