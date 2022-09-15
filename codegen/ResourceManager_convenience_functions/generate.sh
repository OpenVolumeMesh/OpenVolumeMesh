#!/bin/bash

cat vertex_functions.hh > convenience.hh
for Entity in Edge HalfEdge Face HalfFace Cell; 
do
    lower=$(echo ${Entity} | tr '[:upper:]' '[:lower:]')
    sed "s/Vertex/${Entity}/g; s/vertex/${lower}/g" vertex_functions.hh >> convenience.hh
done
