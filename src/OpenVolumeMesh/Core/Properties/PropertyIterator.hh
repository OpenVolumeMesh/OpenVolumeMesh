#pragma once

namespace OpenVolumeMesh {

class PropertyStorageBase;

template<typename Iter>
class PropertyIterator {
public:
    PropertyIterator(Iter _it)
        : it_(_it)
    {}
    PropertyStorageBase * operator*() {
        return it_->get();
    }
    PropertyIterator<Iter> operator++() {
        ++it_;
        return *this;
    }
    bool operator!=(PropertyIterator<Iter> const &other) {
        return it_ != other.it_;
    }
private:
    Iter it_;
};


} // namespace OpenVolumeMesh
