#pragma once
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <type_traits>

namespace OpenVolumeMesh {

template<typename Entity, typename IntT=uint64_t, typename TagT=IntT>
class SmartTagger
{
static_assert(std::is_unsigned<IntT>::value); // overflow must be defined

public:
  using Handle = HandleT<Entity>;
    SmartTagger(ResourceManager const &_mesh, IntT _tag_range = 2)
      : prop_(_mesh.create_private_property<IntT, Entity>("", 0))
      , tag_range_(_tag_range)
    {
      // equality makes the smart tagger pretty useless, but at least it would work
      assert(_tag_range <= std::numeric_limits<IntT>::max());
    }
    void reset() {
      if (current_base_ < std::numeric_limits<IntT>::max() - tag_range_ * 2) {
        current_base_ += tag_range_;
      } else {
        current_base_ = 0;
        prop_.fill(0);
      }
    }
    TagT get(Handle h) const {
      IntT val = prop_[h] - current_base_;
      if (val >= tag_range_) {
        val = 0;
      }
      return static_cast<TagT>(val);
    }

    TagT operator[](Handle h) const {
      return get(h);
    }

    void set(Handle h, TagT tag) {
      prop_[h] = current_base_ + static_cast<IntT>(tag);
    }

    IntT current_base() const {
        return current_base_;
    }


  private:
  // TODO: use private property without shared_ptr indirection!
  PropertyPtr<IntT, Entity> prop_;
  IntT current_base_ = 0;
  IntT tag_range_ = 2; // number of distinct tags
};

template<typename Entity, typename IntT=uint64_t>
using SmartTaggerBool = SmartTagger<Entity, IntT, bool>;

} // namespace OpenVolumeMesh
