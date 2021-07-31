#pragma once

#include <boost/spirit/include/qi.hpp>

#include <boost/bind.hpp>

#include "../MeshGenerator.hpp"

#include <iostream>

namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;

template <typename Iterator>
class tetmesh_grammar : public qi::grammar<Iterator, qi::space_type> {
public:
    explicit tetmesh_grammar(MeshGenerator& _generator) :
        tetmesh_grammar::base_type(content),
        generator_(_generator)
    {}

private:

    MeshGenerator& generator_;

    using Rule = qi::rule<Iterator, qi::space_type>;
    Rule node_section_header {spirit::lit("Vertices") >> qi::int_};
    Rule element_section_header {spirit::lit("Tetrahedra") >>
                qi::int_[boost::bind(&MeshGenerator::set_num_cells, &generator_, ::_1)]};
    Rule node {qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               qi::double_};
    Rule element {qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)]};
    Rule content {node_section_header >> *node >> element_section_header >> *element};
};
