#pragma once
#include "verification/Sha256.h"
#include <ostream>
#include <streambuf>
#include <stdexcept>

namespace verification {
// Uses the repository SHA-256 unchanged. Only 8 KiB of formatted bytes are
// buffered; digest() snapshots the hash without discarding or retaining input.
class Sha256Stream {
    class Buffer final : public std::streambuf {
        Sha256 hash_;
        std::array<char,8192> bytes_{};
        void drain() {
            hash_.update(reinterpret_cast<const unsigned char*>(pbase()),
                         static_cast<std::size_t>(pptr()-pbase()));
            setp(bytes_.data(),bytes_.data()+bytes_.size());
        }
    protected:
        int_type overflow(int_type ch) override {
            drain();
            if(!traits_type::eq_int_type(ch,traits_type::eof())) {
                *pptr()=traits_type::to_char_type(ch);pbump(1);
            }
            return traits_type::not_eof(ch);
        }
        int sync() override {drain();return 0;}
    public:
        Buffer(){setp(bytes_.data(),bytes_.data()+bytes_.size());}
        std::string digest() {drain();auto snapshot=hash_;return snapshot.finish();}
    } buffer_;
    std::ostream stream_{&buffer_};
public:
    Sha256Stream()=default;
    Sha256Stream(const Sha256Stream&)=delete;
    Sha256Stream& operator=(const Sha256Stream&)=delete;
    std::ostream& stream() noexcept {return stream_;}
    std::string digest() {
        if(!stream_ || stream_.rdbuf()!=&buffer_)
            throw std::runtime_error("invalid SHA256 output stream");
        stream_.flush();
        if(!stream_)throw std::runtime_error("SHA256 output stream flush failed");
        return buffer_.digest();
    }
};
}
