#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <optional>

namespace OpenVolumeMesh {

class TriangleTopology;

/**
 * Consistent labeling of the entities a tet is comprised of.
 * ABC is an interior half-face (part of the cell,
 * CCW oriented as seen from "inside" the cell),
 * with D as its opposite vertex.
 */
class OVM_EXPORT TetTopology {
public:
    enum VertexLabel : uint8_t {
        A, B, C, D
    };
    enum HalfEdgeLabel : uint8_t {
        // Note: halfedge i+3 is vertex-disjunct with halfedge i (\in [0..2])
        AB, BC, CA,
        CD, AD, BD,
        // Opposites, with flip of bit 3:
        BA = 8 + AB, CB = 8 + BC, AC = 8 + CA,
        DC = 8 + CD, DA = 8 + AD, DB = 8 + BD,
    };

    static constexpr HalfEdgeLabel opposite(HalfEdgeLabel _hel) {
        return static_cast<HalfEdgeLabel>(_hel ^ 8);
    }

    /// Is a halfedge label 'forward', i.e. an index into heh_?
    static constexpr bool is_forward(HalfEdgeLabel _hel) {
        return (_hel & 8) == 0;
    }

    template<VertexLabel From, VertexLabel To>
    static constexpr HalfEdgeLabel hel() {
        static_assert(From != To);
        // regular
        if constexpr(From == A && To == B) { return AB; }
        if constexpr(From == B && To == C) { return BC; }
        if constexpr(From == C && To == A) { return CA; }
        if constexpr(From == C && To == D) { return CD; }
        if constexpr(From == A && To == D) { return AD; }
        if constexpr(From == B && To == D) { return BD; }
        // opposite
        if constexpr(From == B && To == A) { return BA; }
        if constexpr(From == C && To == B) { return CB; }
        if constexpr(From == A && To == C) { return AC; }
        if constexpr(From == D && To == C) { return DC; }
        if constexpr(From == D && To == A) { return DA; }
        if constexpr(From == D && To == B) { return DB; }
    }

    template<HalfEdgeLabel HEL>
    static constexpr VertexLabel hel_from() {
        if constexpr(HEL == AB || HEL == AC || HEL == AD) { return A; }
        if constexpr(HEL == BA || HEL == BC || HEL == BD) { return B; }
        if constexpr(HEL == CA || HEL == CB || HEL == CD) { return C; }
        if constexpr(HEL == DA || HEL == DB || HEL == DC) { return D; }
    }

    template<HalfEdgeLabel HEL>
    static constexpr VertexLabel hel_to() {
        if constexpr(HEL == BA || HEL == CA || HEL == DA) { return A; }
        if constexpr(HEL == AB || HEL == CB || HEL == DB) { return B; }
        if constexpr(HEL == AC || HEL == BC || HEL == DC) { return C; }
        if constexpr(HEL == AD || HEL == BD || HEL == CD) { return D; }
    }

    /// Bits [2..3] indicate which HF, Bits [0..1] indicate start vertex
    enum HalfFaceLabel : uint8_t {
        OppA =  0, BDC = OppA+1, CBD = OppA+2, DCB = OppA+3,
        OppB =  4, ACD = OppB+1, CDA = OppB+2, DAC = OppB+3,
        OppC =  8, ADB = OppC+1, BAD = OppC+2, DBA = OppC+3,
        OppD = 12, ABC = OppD+1, BCA = OppD+2, CAB = OppD+3,

    /// Values are HalfFaceLabel of the opposite HF with bit 4 set
        OuterOppA = OppA+16, BCD = BDC+16, CDB = CBD+16, DBC = DCB+16,
        OuterOppB = OppB+16, ADC = ACD+16, CAD = CDA+16, DCA = DAC+16,
        OuterOppC = OppC+16, ABD = ADB+16, BDA = BAD+16, DAB = DBA+16,
        OuterOppD = OppD+16, ACB = ABC+16, BAC = BCA+16, CBA = CAB+16,
    };

    template<HalfFaceLabel HFL>
    static constexpr HalfFaceLabel opposite() {
        return static_cast<HalfFaceLabel>(HFL ^ 16);
    }

    static constexpr HalfFaceLabel opposite(HalfFaceLabel _hel) {
        return static_cast<HalfFaceLabel>(_hel ^ 16);
    }

    static constexpr bool is_inner(HalfFaceLabel _hfl) {
        return (_hfl & 16) == 0;
    }

    template<HalfFaceLabel HFL>
    static constexpr HalfFaceLabel inner() {
        return static_cast<HalfFaceLabel>(HFL & ~16);
    }
    static constexpr HalfFaceLabel inner(HalfFaceLabel _hel) {
        return static_cast<HalfFaceLabel>(_hel & ~16);
    }

    template<HalfFaceLabel HFL>
    static constexpr HalfFaceLabel outer() {
        return static_cast<HalfFaceLabel>(HFL | 16);
    }
    static constexpr HalfFaceLabel outer(HalfFaceLabel _hel) {
        return static_cast<HalfFaceLabel>(_hel | 16);
    }
    static constexpr bool has_start(HalfFaceLabel _hfl) {
        return _hfl != OppA
            && _hfl != OppB
            && _hfl != OppC
            && _hfl != OppD
            && _hfl != OuterOppA
            && _hfl != OuterOppB
            && _hfl != OuterOppC
            && _hfl != OuterOppD;
    }

    /// Get the label of the Idx'th vertex in the HalfFace HFL
    template<HalfFaceLabel HFL, uint8_t Idx>
    static constexpr VertexLabel hfl_vl() {
        static_assert(Idx < 3);
        // OppX have no defined starting vertex
        static_assert(has_start(HFL));
        if constexpr (Idx == 0) {
            if constexpr (HFL == ACD || HFL == ADB || HFL == ABC) {return A;}
            if constexpr (HFL == ADC || HFL == ABD || HFL == ACB) {return A;}
            if constexpr (HFL == BDC || HFL == BAD || HFL == BCA) {return B;}
            if constexpr (HFL == BCD || HFL == BDA || HFL == BAC) {return B;}
            if constexpr (HFL == CBD || HFL == CDA || HFL == CAB) {return C;}
            if constexpr (HFL == CDB || HFL == CAD || HFL == CBA) {return C;}
            if constexpr (HFL == DCB || HFL == DAC || HFL == DBA) {return D;}
            if constexpr (HFL == DBC || HFL == DCA || HFL == DAB) {return D;}
        } else if constexpr (Idx == 1) {
            if constexpr (HFL == DAC || HFL == BAD || HFL == CAB) {return A;}
            if constexpr (HFL == CAD || HFL == DAB || HFL == BAC) {return A;}
            if constexpr (HFL == CBD || HFL == DBA || HFL == ABC) {return B;}
            if constexpr (HFL == DBC || HFL == ABD || HFL == CBA) {return B;}
            if constexpr (HFL == DCB || HFL == ACD || HFL == BCA) {return C;}
            if constexpr (HFL == BCD || HFL == DCA || HFL == ACB) {return C;}
            if constexpr (HFL == BDC || HFL == CDA || HFL == ADB) {return D;}
            if constexpr (HFL == CDB || HFL == ADC || HFL == BDA) {return D;}
        } else if constexpr (Idx == 2) {
            if constexpr (HFL == CDA || HFL == DBA || HFL == BCA) {return A;}
            if constexpr (HFL == DCA || HFL == BDA || HFL == CBA) {return A;}
            if constexpr (HFL == DCB || HFL == ADB || HFL == CAB) {return B;}
            if constexpr (HFL == CDB || HFL == DAB || HFL == ACB) {return B;}
            if constexpr (HFL == BDC || HFL == DAC || HFL == ABC) {return C;}
            if constexpr (HFL == DBC || HFL == ADC || HFL == BAC) {return C;}
            if constexpr (HFL == CBD || HFL == ACD || HFL == BAD) {return D;}
            if constexpr (HFL == BCD || HFL == CAD || HFL == ABD) {return D;}
        }
        return static_cast<VertexLabel>(-1); // unreachable
    }
    /// Get the label of the Idx'th halfedge in the HalfFace HFL
    template<HalfFaceLabel HFL, uint8_t Idx>
    static constexpr HalfEdgeLabel hfl_hel() {
        constexpr VertexLabel from = hfl_vl<HFL, Idx>();
        constexpr VertexLabel to   = hfl_vl<HFL, (Idx+1)%3>();
        return hel<from,to>();
    }

    /// Find a vertex labeling for the cell `ch` with the given
    /// halfface abc and vertex a.
    ///
    /// `abc` must be a a halfface of a the topological tetrahedron `ch`,
    /// and contain the vertex `a`.
    TetTopology(TopologyKernel const &mesh, CH ch, HFH abc, VH a=VH());

    /// Find a vertex labeling for the cell incident to `abc`, such that
    /// `a` is the given vertex handle, and a, b, and c form `abc`.
    ///
    /// `ch` MUST be a topological tetrahedron
    /// `hfh` MUST be a halface of a topological tetrahedron and contain `a`.
    TetTopology(TopologyKernel const &mesh, HFH abc, VH a=VH());

    /// Find a vertex labeling for the cell `ch`, such that `a` is the given
    /// vertex.
    /// `ch` MUST be a topological tetrahedron
    /// `a` MUST be a vertex of `ch`.
    TetTopology(TopologyKernel const &mesh, CH ch, VH a);

    /// Find any vertex labeling for the tetrahedral `ch`.
    /// `ch` MUST be a topological tetrahedron.
    TetTopology(TopologyKernel const &mesh, CH ch);

    bool operator==(TetTopology const& other) const;

    const auto &halfface_handles() const {return hfh_;}

    template<VertexLabel I>
    constexpr VH vh() const {return vh_[I];}
    //constexpr VH vh(VertexLabel i) const {return vh_[i];}

    template<HalfEdgeLabel I>
    constexpr HEH heh() const {
        if constexpr(is_forward(I)) {
            return heh_[I];
        } else {
            return heh_[I & 7].opposite_handle();
        }
    }
    //constexpr HEH heh(HalfEdgeLabel i) const {return heh_[i];}

    template<HalfFaceLabel I>
    constexpr HFH hfh() const {
        if constexpr(is_inner(I)) {
            return hfh_[I>>2];
        } else {
            return hfh_[(I & 15)>>2].opposite_handle();
        }
    }
    //constexpr HFH hfh(HalfFaceLabel i) const {return hfh_[i>>2];}

    template<VertexLabel X, VertexLabel Y>
    constexpr HEH heh() const {
        return heh<hel<X,Y>()>();
    }

    // convenience helper functions:

    constexpr VH a() const {return vh<A>();}
    constexpr VH b() const {return vh<B>();}
    constexpr VH c() const {return vh<C>();}
    constexpr VH d() const {return vh<D>();}

    // half-edges
    constexpr HEH ab() const {return heh<AB>();}
    constexpr HEH ca() const {return heh<CA>();}
    constexpr HEH ad() const {return heh<AD>();}
    constexpr HEH bc() const {return heh<BC>();}
    constexpr HEH bd() const {return heh<BD>();}
    constexpr HEH cd() const {return heh<CD>();}

    constexpr HEH ac() const {return ca().opposite_handle();}
    constexpr HEH ba() const {return ab().opposite_handle();}
    constexpr HEH da() const {return ad().opposite_handle();}
    constexpr HEH dc() const {return cd().opposite_handle();}
    constexpr HEH db() const {return bd().opposite_handle();}
    constexpr HEH cb() const {return bc().opposite_handle();}

    // inner half-faces with all permutations
    constexpr HFH bdc() const {return hfh<BDC>();}
    constexpr HFH dcb() const {return hfh<DCB>();}
    constexpr HFH cbd() const {return hfh<CBD>();}

    constexpr HFH acd() const {return hfh<ACD>();}
    constexpr HFH cda() const {return hfh<CDA>();}
    constexpr HFH dac() const {return hfh<DAC>();}

    constexpr HFH bad() const {return hfh<BAD>();}
    constexpr HFH adb() const {return hfh<ADB>();}
    constexpr HFH dba() const {return hfh<DBA>();}

    constexpr HFH abc() const {return hfh<ABC>();}
    constexpr HFH bca() const {return hfh<BCA>();}
    constexpr HFH cab() const {return hfh<CAB>();}
    // outer hfhs with all permutations
    constexpr HFH bcd() const {return hfh<BCD>();}
    constexpr HFH dbc() const {return hfh<DBC>();}
    constexpr HFH cdb() const {return hfh<CDB>();}

    constexpr HFH adc() const {return hfh<ADC>();}
    constexpr HFH cad() const {return hfh<CAD>();}
    constexpr HFH dca() const {return hfh<DCA>();}

    constexpr HFH bda() const {return hfh<BDA>();}
    constexpr HFH abd() const {return hfh<ABD>();}
    constexpr HFH dab() const {return hfh<DAB>();}

    constexpr HFH acb() const {return hfh<ACB>();}
    constexpr HFH bac() const {return hfh<BAC>();}
    constexpr HFH cba() const {return hfh<CBA>();}

    /** Create a TriangleTopology for the given halfface.
     * triangle_topology<XYZ>() return a TriangleTopology,
     * such that a=X, b=Y, c=Z.
     */
    template<HalfFaceLabel I>
    constexpr TriangleTopology triangle_topology() const;

    /** Create a TriangleTopology for the given inner halfface.
     * Prefer the templated version triangle_topology<HalfFaceLabel> when
     * the half-face label is known at compile time.
     */
    TriangleTopology triangle_topology(HalfFaceLabel i) const;
    
    // todo:
    //  - rotated<HalfFaceLabel>, rotated(HalfFaceLabel)
    //  - rotated(hfh), rotated(hfh, vh)


    std::optional<VertexLabel> get_label(VH _vh) const;
    std::optional<HalfEdgeLabel> get_label(HEH _heh) const;
    std::optional<HalfFaceLabel> get_label(HFH _hfh) const;
    std::optional<HalfFaceLabel> get_label(HFH _hfh, VH _first) const;

private:
    template<HalfFaceLabel I>
    constexpr HFH& hfh() {static_assert(I<16); return hfh_[I>>2];}

    std::array<VH, 4> vh_;
    std::array<HEH, 6> heh_;
    std::array<HFH, 4> hfh_;
};


} // namespace OpenVolumeMesh

// TODO: implicit instantiations to improve compile times.
#include <OpenVolumeMesh/Unstable/Topology/TetTopology_impl.hh>
