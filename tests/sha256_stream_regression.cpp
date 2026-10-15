#if __has_include("verification/Sha256Stream.h")
#include "verification/Sha256Stream.h"
#include <iostream>
#include <stdexcept>
#include <limits>

void require(bool ok){if(!ok)throw std::runtime_error("SHA256 stream regression");}
std::string plain(const std::string& s){verification::Sha256 h;h.update(reinterpret_cast<const unsigned char*>(s.data()),s.size());return h.finish();}
int main(){try{
    verification::Sha256Stream empty;
    require(empty.digest()=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    verification::Sha256Stream h;auto& o=h.stream();o<<"abc";
    require(h.digest()=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    require(h.digest()==plain("abc"));o<<'d';require(h.digest()==plain("abcd"));
    verification::Sha256Stream million;
    for(int i=0;i<1000000;++i)million.stream().put('a');
    require(million.digest()=="cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
    for(int n:{55,56,63,64,65,8191,8192,8193,32771}){
        std::string s(n,'x');s[n/2]='\0';verification::Sha256Stream x;
        x.stream().write(s.data(),s.size());x.stream()<<std::flush;require(x.digest()==plain(s));
    }
    verification::Sha256Stream formatted;std::ostringstream reference;
    auto write=[](std::ostream& x){x.imbue(std::locale::classic());x<<std::hexfloat<<-0.<<','<<1./3<<','<<std::numeric_limits<double>::infinity()<<','<<std::quoted("escape\"\\")<<std::setw(12)<<17;};
    write(formatted.stream());write(reference);require(formatted.digest()==plain(reference.str()));
    formatted.stream().setstate(std::ios::badbit);
    bool rejected=false;try{formatted.digest();}catch(const std::exception&){rejected=true;}require(rejected);
    std::cout<<"SHA256 stream known vectors, boundaries, embedded NUL, formatting and failed stream passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#else
#include <iostream>
int main(){std::cerr<<"missing bounded SHA256 stream sink\n";return 1;}
#endif
