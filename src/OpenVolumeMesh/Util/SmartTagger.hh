#pragma once
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>

namespace OpenVolumeMesh {

template<typename Entity, typename IntT=uint64_t>
class SmartTagger
{
  using Handle = HandleT<Entity>;
  public:
    SmartTagger(ResourceManager const &_mesh)
      : prop_(_mesh.create_private_property<IntT, Entity>("", 0))
      , generation_(1)
    {}
    void reset() {
      ++generation_;
      if (generation_ == 0) {
        prop_.fill(0);
        ++generation_;
      }
    }
    void tag(Handle h) {
      prop_[h] = generation_;
    }

    void untag(Handle h) {
      prop_[h] = 0;
    }

    bool operator[](Handle h) const {
      return prop_[h] == generation_;
    }

  using T = IntT;
  private:
  // TODO: use private property without shared_ptr indirection!
  PropertyPtr<IntT, Entity> prop_;
  T generation_;
};

} // namespace OpenVolumeMesh
