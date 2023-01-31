#include <OpenVolumeMesh/Unstable/Topology/TetTopology.hh>
#include <OpenVolumeMesh/Unstable/Topology/TriangleTopology.hh>
#include <type_traits>

namespace OpenVolumeMesh {

static HFH find_halfface(TopologyKernel const &_mesh, CH _ch, VH _vh) {
    for (const auto hfh: _mesh.cell_halffaces(_ch)) {
        for (const auto vh: _mesh.halfface_vertices(hfh)) {
            if (vh == _vh)
                return hfh;
        }
    }
    assert(false);
    return {};
}


TetTopology::TetTopology(const TopologyKernel &mesh,
                         CH _ch,
                         HFH _abc,
                         VH _a)
{
    auto abc_he_it = mesh.hfhe_iter(_abc, 2);
    if (_a.is_valid()) {
        while(mesh.from_vertex_handle(*abc_he_it) != _a) {
            ++abc_he_it;
        }
    }
    heh_[AB] = *abc_he_it++;
    heh_[BC] = *abc_he_it++;
    heh_[CA] = *abc_he_it;

    vh_[A] = mesh.from_vertex_handle(ab());
    vh_[B] = mesh.from_vertex_handle(bc());
    vh_[C] = mesh.from_vertex_handle(ca());

    assert(b() == mesh.to_vertex_handle(ab()));
    assert(c() == mesh.to_vertex_handle(bc()));
    assert(a() == mesh.to_vertex_handle(ca()));

    hfh<ABC>() = _abc;
    for (const auto cur_hfh: mesh.cell_halffaces(_ch)) {
        if (cur_hfh == _abc) continue;
        auto heh_it = mesh.hfhe_iter(cur_hfh, 2);
        while (heh_it) {
            auto heh = *(heh_it++);
            if (heh == ba()) {
                hfh<BAD>() = cur_hfh;
                heh_[AD] = *heh_it;
                break;
            }
            if (heh == cb()) {
                hfh<CBD>() = cur_hfh;
                heh_[BD] = *heh_it;
                break;
            }
            if (heh == ac()) {
                hfh<ACD>() = cur_hfh;
                heh_[CD] = *heh_it;
                break;
            }
        }
    }
    vh_[D] = mesh.to_vertex_handle(ad());

    assert(d() != a());
    assert(d() != b());
    assert(d() != c());

    assert(d() == mesh.to_vertex_handle(ad()));
    assert(d() == mesh.to_vertex_handle(bd()));
    assert(d() == mesh.to_vertex_handle(cd()));

    assert(a() == mesh.from_vertex_handle(ad()));
    assert(b() == mesh.from_vertex_handle(bd()));
    assert(c() == mesh.from_vertex_handle(cd()));
}


TetTopology::TetTopology(TopologyKernel const &_mesh,
        HFH _abc, VH _a)
    : TetTopology(_mesh, _mesh.incident_cell(_abc), _abc, _a)
{}

TetTopology::TetTopology(const TopologyKernel &_mesh,
                         CH _ch, VH _a)
    : TetTopology(_mesh, _ch, find_halfface(_mesh, _ch, _a), _a)
{}

TetTopology::TetTopology(const TopologyKernel &_mesh, CH _ch)
    : TetTopology(_mesh, _ch, *_mesh.chf_iter(_ch), VH())
{}

bool TetTopology::operator==(TetTopology const& other) const
{
  return vh_ == other.vh_
    && heh_ == other.heh_
    && hfh_ == other.hfh_;
}

auto TetTopology::get_label(VH _vh) const
    -> std::optional<VertexLabel> 
{
    for (uint8_t idx = 0; idx < vh_.size(); ++idx) {
        if (_vh == vh_[idx]) {
            return static_cast<VertexLabel>(idx);
        }
    }
    return {};
}
    auto TetTopology::get_label(HEH _heh) const
    -> std::optional<HalfEdgeLabel>
{
    const auto eh = _heh.edge_handle();
    for (uint8_t idx = 0; idx < heh_.size(); ++idx) {
        if (heh_[idx].edge_handle() == eh) {
            auto hel = static_cast<HalfEdgeLabel>(idx);
            if (heh_[idx] == _heh) {
                return hel;
            } else {
                return opposite(hel);
            }
        }
    }
    return {};
}
auto TetTopology::get_label(HFH _hfh) const
    -> std::optional<HalfFaceLabel>
{
    const auto fh = _hfh.face_handle();
    for (uint8_t idx = 0; idx < hfh_.size(); ++idx) {
        if (hfh_[idx].face_handle() == fh) {
            auto hfl = static_cast<HalfFaceLabel>(idx << 2);
            if (hfh_[idx] == _hfh) {
                return hfl;
            } else  {
                return opposite(hfl);
            }
        }
    }
    return {};
}
namespace detail {
    template<TetTopology::HalfFaceLabel HFL>
    std::optional<TetTopology::HalfFaceLabel> try_get_label(TetTopology const& tettopo, HFH _hfh, TetTopology::VertexLabel _first)
    {
        using HalfFaceLabel = TetTopology::HalfFaceLabel;
        if (tettopo.hfh<HFL>() == _hfh) {
            if (tettopo.hfl_vl<static_cast<HalfFaceLabel>(HFL+1), 0>() == _first) { return static_cast<HalfFaceLabel>(HFL+1); }
            if (tettopo.hfl_vl<static_cast<HalfFaceLabel>(HFL+2), 0>() == _first) { return static_cast<HalfFaceLabel>(HFL+2); }
            if (tettopo.hfl_vl<static_cast<HalfFaceLabel>(HFL+3), 0>() == _first) { return static_cast<HalfFaceLabel>(HFL+3); }
        } else if (tettopo.hfh<HFL>() == _hfh.opposite_handle()) {
            return try_get_label<TetTopology::opposite<HFL>()>(tettopo, _hfh, _first);
        }
        return {};
    }
}
auto TetTopology::get_label(HFH _hfh, VH _first) const
    -> std::optional<HalfFaceLabel>
{
    auto first_label = get_label(_first);
    if (!first_label)
        return {};
    if (auto lbl = detail::try_get_label<OppA>(*this, _hfh, *first_label)) { return *lbl;}
    if (auto lbl = detail::try_get_label<OppB>(*this, _hfh, *first_label)) { return *lbl;}
    if (auto lbl = detail::try_get_label<OppC>(*this, _hfh, *first_label)) { return *lbl;}
    if (auto lbl = detail::try_get_label<OppD>(*this, _hfh, *first_label)) { return *lbl;}
    return {};
}
namespace detail {
    template<TetTopology::HalfFaceLabel HFL>
    using hfl_dummy = std::integral_constant<TetTopology::HalfFaceLabel, HFL>;
    template<typename F>
    constexpr auto dynamic_dispatch(F fun, TetTopology::HalfFaceLabel hfl) {
        switch (hfl) {
            case TetTopology::OppA:      return fun(hfl_dummy<TetTopology::OppA>());
            case TetTopology::BDC:       return fun(hfl_dummy<TetTopology::BDC>());
            case TetTopology::CBD:       return fun(hfl_dummy<TetTopology::CBD>());
            case TetTopology::DCB:       return fun(hfl_dummy<TetTopology::DCB>());
            case TetTopology::OppB:      return fun(hfl_dummy<TetTopology::OppB>());
            case TetTopology::ACD:       return fun(hfl_dummy<TetTopology::ACD>());
            case TetTopology::CDA:       return fun(hfl_dummy<TetTopology::CDA>());
            case TetTopology::DAC:       return fun(hfl_dummy<TetTopology::DAC>());
            case TetTopology::OppC:      return fun(hfl_dummy<TetTopology::OppC>());
            case TetTopology::ADB:       return fun(hfl_dummy<TetTopology::ADB>());
            case TetTopology::BAD:       return fun(hfl_dummy<TetTopology::BAD>());
            case TetTopology::DBA:       return fun(hfl_dummy<TetTopology::DBA>());
            case TetTopology::OppD:      return fun(hfl_dummy<TetTopology::OppD>());
            case TetTopology::ABC:       return fun(hfl_dummy<TetTopology::ABC>());
            case TetTopology::BCA:       return fun(hfl_dummy<TetTopology::BCA>());
            case TetTopology::CAB:       return fun(hfl_dummy<TetTopology::CAB>());
            case TetTopology::OuterOppA: return fun(hfl_dummy<TetTopology::OuterOppA>());
            case TetTopology::BCD:       return fun(hfl_dummy<TetTopology::BCD>());
            case TetTopology::CDB:       return fun(hfl_dummy<TetTopology::CDB>());
            case TetTopology::DBC:       return fun(hfl_dummy<TetTopology::DBC>());
            case TetTopology::OuterOppB: return fun(hfl_dummy<TetTopology::OuterOppB>());
            case TetTopology::ADC:       return fun(hfl_dummy<TetTopology::ADC>());
            case TetTopology::CAD:       return fun(hfl_dummy<TetTopology::CAD>());
            case TetTopology::DCA:       return fun(hfl_dummy<TetTopology::DCA>());
            case TetTopology::OuterOppC: return fun(hfl_dummy<TetTopology::OuterOppC>());
            case TetTopology::ABD:       return fun(hfl_dummy<TetTopology::ABD>());
            case TetTopology::BDA:       return fun(hfl_dummy<TetTopology::BDA>());
            case TetTopology::DAB:       return fun(hfl_dummy<TetTopology::DAB>());
            case TetTopology::OuterOppD: return fun(hfl_dummy<TetTopology::OuterOppD>());
            case TetTopology::ACB:       return fun(hfl_dummy<TetTopology::ACB>());
            case TetTopology::BAC:       return fun(hfl_dummy<TetTopology::BAC>());
            case TetTopology::CBA:       return fun(hfl_dummy<TetTopology::CBA>());
        }
        throw std::runtime_error("invalid halfface label");
    }
}
TriangleTopology TetTopology::triangle_topology(HalfFaceLabel _hfl) const {
    return detail::dynamic_dispatch([this](auto dummy){
            constexpr auto hfl = decltype(dummy)::value;
            if constexpr (has_start(hfl)) {
                return triangle_topology<hfl>();
            } else {
                throw std::runtime_error("invalid hfl");
                // HACK: never reachable, make compiler happy by returning right type:
                return triangle_topology<ABC>();
            }
            }, _hfl);
}


} // namespace OpenVolumeMesh
