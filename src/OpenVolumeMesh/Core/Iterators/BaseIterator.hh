#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>

#include <iterator>

namespace OpenVolumeMesh {
class TopologyKernel;

template <
class OH /* Output handle type */>
class OVM_EXPORT BaseIterator {
public:

    // STL compliance
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef int						        difference_type;
    typedef const OH				        value_type;
    typedef const OH*				        pointer;
    typedef const OH&				        reference;


    BaseIterator(const TopologyKernel* _mesh, const OH& _ch) :
        valid_(true), cur_handle_(_ch), mesh_(_mesh) {}

    explicit BaseIterator(const TopologyKernel* _mesh) :
        valid_(true), mesh_(_mesh) {}

    // STL compliance (needs to have default constructor)
    BaseIterator() : valid_(false), mesh_(nullptr) {}
    BaseIterator(const BaseIterator& _c) = default;
    virtual ~BaseIterator() = default;

    BaseIterator& operator=(const BaseIterator& _c) = default;

    bool operator== (const BaseIterator& _c) const {
        return (this->cur_handle_ == _c.cur_handle() &&
                this->valid_ == _c.valid() &&
                this->mesh_ == _c.mesh());
    }
    bool operator!= (const BaseIterator& _c) const {
        return !this->operator==(_c);
    }

    pointer operator->() const {
        return &cur_handle_;
    }

    reference operator*() const {
        return cur_handle_;
    }

    bool operator< (const BaseIterator& _c) const {
        return cur_handle_.idx() < _c.cur_handle_.idx();
    }

    operator bool() const {
        return valid_;
    }

    void valid(bool _valid) {
        valid_ = _valid;
    }
    // TODO: [[deprecated("Use is_valid() instead")]]
    bool valid() const {
        return valid_;
    }
    bool is_valid() const {
        return valid_;
    }
    void cur_handle(const OH& _h) {
        cur_handle_ = _h;
    }
    reference cur_handle() const {
        return cur_handle_;
    }
    const TopologyKernel* mesh() const {
        return mesh_;
    }

private:

    bool valid_;
    OH cur_handle_;
    const TopologyKernel* mesh_;
};

template<class I>
using is_ovm_iterator = std::is_base_of<BaseIterator<typename std::remove_const<typename I::value_type>::type>, I>;

} // namespace OpenVolumeMesh
