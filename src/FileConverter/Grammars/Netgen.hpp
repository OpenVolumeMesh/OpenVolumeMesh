#pragma once

#include <boost/spirit/include/qi.hpp>

#include <boost/bind.hpp>

#include "../MeshGenerator.hpp"

#include <iostream>

namespace qi = boost::spirit::qi;
namespace spirit = boost::spirit;

void print() {
    static int c = 0;
    std::cerr << "Hey there!" << c++ << std::endl;
}

template <typename Iterator>
class netgen_grammar : public qi::grammar<Iterator/*, qi::space_type*/> {
public:
    explicit netgen_grammar(MeshGenerator& _generator) :
        netgen_grammar::base_type(content),
        generator_(_generator)
    {}

private:
    MeshGenerator& generator_;

    using Rule = qi::rule<Iterator>;

    Rule space {spirit::ascii::space - spirit::eol};
    Rule node_section_header {qi::int_ /* Number of vertices */ >> spirit::eol};
    Rule element_section_header {qi::int_[boost::bind(&MeshGenerator::set_num_cells, &generator_, ::_1)] /* Number of tetrahedra */ >> spirit::eol};
    Rule face_section_header {qi::int_ /* Number of faces */ >> spirit::eol};
    Rule node {*space >> qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               *space >> qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >>
               *space >> qi::double_[boost::bind(&MeshGenerator::add_vertex_component, &generator_, ::_1)] >> spirit::eol};
    Rule element {*space >> qi::int_ >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  *space >> qi::int_[boost::bind(&MeshGenerator::add_cell_vertex, &generator_, ::_1)] >>
                  spirit::eol};
    Rule face {*space >> qi::int_ >>
               *space >> qi::int_ >>
               *space >> qi::int_ >>
               *space >> qi::int_ >>
               spirit::eol};
    Rule content {node_section_header >> *node >>
                  element_section_header >> *element >>
                  face_section_header >> *face};

};
