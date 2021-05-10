#pragma once

#ifdef __OSX__
# include <OpenGL/gl3.h>
#else
# include <GL3/gl3.h>
#endif
#include <morph/tools.h>
#include <morph/VisualModel.h>
#include <morph/ColourMap.h>
#include <morph/Scale.h>
#include <morph/MathAlgo.h>
#include <morph/Vector.h>
#include <iostream>
#include <vector>
#include <array>

#include "tissue.h"

template <class T>
class retvisual : public morph::VisualModel
{
public:
    //! Single constructor for simplicity
    retvisual(GLuint sp, GLuint tsp, const retina<T>* _r, const morph::Vector<float> _offset)
    {
        this->shaderprog = sp;
        this->tshaderprog = tsp;
        this->mv_offset = _offset;
        this->viewmatrix.translate (this->mv_offset);
        // Defaults for z and colourScale
        this->zScale.setParams (1, 0);
        this->colourScale.do_autoscale = true;
        this->ret = _r;
        // Note: VisualModel::finalize() should be called before rendering
    }

    //! Initialize as a rectangle made of 4 triangles for each rect, with z position
    //! of each of the 4 outer edges of the triangles interpolated, but a single colour
    //! for each rectangle. Gives a smooth surface in which you can see the pixels.
    void initializeVertices()
    {
        float dx = (float)this->ret->dx[0];
        float hx = 0.5f * dx;
        float dy = (float)this->ret->dx[1];
        float vy = 0.5f * dy;

        unsigned int idx = 0;

        // This changes to use the interaction parameter of the retina object.
        this->dcopy.resize (this->ret->interaction.size());
        this->dcolour.resize (this->ret->interaction.size());
        this->dcolour2.resize (this->ret->interaction.size());
        for (unsigned int i = 0; i < this->ret->interaction.size(); ++i) {
            std::cout << "ret->interaction["<<i<<"] = " << this->ret->interaction[i] << std::endl;
            this->dcolour[i] = this->ret->interaction[i][0];
            this->dcolour2[i] = this->ret->interaction[i][1];
        }
        this->colourScale.transform (this->dcolour, this->dcolour);
        std::pair<T,T> maxmin = morph::MathAlgo::maxmin (this->dcolour);
        this->colourScale.autoscaled = false;
        this->colourScale.transform (this->dcolour2, this->dcolour2);
        std::pair<T,T> maxmin2 = morph::MathAlgo::maxmin (this->dcolour2);
        std::cout << "R maxmin: " << maxmin.first << ","<< maxmin.second
                  << " and G maxmin: " << maxmin2.first << ","<< maxmin2.second << std::endl;

        float z = 0.0f;

        morph::Vector<float> vtx_0, vtx_1, vtx_2;
        for (unsigned int ri = 0; ri < this->ret->num(); ++ri) {

            std::array<float, 3> clr = this->setColour (ri);

            // First push the 5 positions of the triangle vertices, starting with the centre. All at z=0
            this->vertex_push (this->ret->posn[ri][0], this->ret->posn[ri][1], z, this->vertexPositions);
            // Use the centre position as the first location for finding the normal vector
            vtx_0 = {{this->ret->posn[ri][0], this->ret->posn[ri][1], z}};
            // NE vertex
            this->vertex_push (this->ret->posn[ri][0]+hx, this->ret->posn[ri][1]+vy, z, this->vertexPositions);
            vtx_1 = {{this->ret->posn[ri][0]+hx, this->ret->posn[ri][1]+vy, z}};
            // SE vertex
            this->vertex_push (this->ret->posn[ri][0]+hx, this->ret->posn[ri][1]-vy, z, this->vertexPositions);
            vtx_2 = {{this->ret->posn[ri][0]+hx, this->ret->posn[ri][1]-vy, z}};
            // SW vertex
            this->vertex_push (this->ret->posn[ri][0]-hx, this->ret->posn[ri][1]-vy, z, this->vertexPositions);
            // NW vertex
            this->vertex_push (this->ret->posn[ri][0]-hx, this->ret->posn[ri][1]+vy, z, this->vertexPositions);

            // From vtx_0,1,2 compute normal. This sets the correct normal, but note
            // that there is only one 'layer' of vertices; the back of the
            // HexGridVisual will be coloured the same as the front. To get lighting
            // effects to look really good, the back of the surface could need the
            // opposite normal.
            morph::Vector<float> plane1 = vtx_1 - vtx_0;
            morph::Vector<float> plane2 = vtx_2 - vtx_0;
            morph::Vector<float> vnorm = plane2.cross (plane1);
            vnorm.renormalize();
            this->vertex_push (vnorm, this->vertexNormals);
            this->vertex_push (vnorm, this->vertexNormals);
            this->vertex_push (vnorm, this->vertexNormals);
            this->vertex_push (vnorm, this->vertexNormals);
            this->vertex_push (vnorm, this->vertexNormals);

            // Five vertices with the same colour
            this->vertex_push (clr, this->vertexColors);
            this->vertex_push (clr, this->vertexColors);
            this->vertex_push (clr, this->vertexColors);
            this->vertex_push (clr, this->vertexColors);
            this->vertex_push (clr, this->vertexColors);

            // Define indices now to produce the 4 triangles in the hex
            this->indices.push_back (idx+1);
            this->indices.push_back (idx);
            this->indices.push_back (idx+2);

            this->indices.push_back (idx+2);
            this->indices.push_back (idx);
            this->indices.push_back (idx+3);

            this->indices.push_back (idx+3);
            this->indices.push_back (idx);
            this->indices.push_back (idx+4);

            this->indices.push_back (idx+4);
            this->indices.push_back (idx);
            this->indices.push_back (idx+1);

            idx += 5; // 5 vertices (each of 3 floats for x/y/z), 15 indices.
        }
    }

    morph::ColourMap<float> cm;
    morph::Scale<T, float> colourScale;
    morph::Scale<T, float> colourScale2;
    morph::Scale<T, float> zScale;

protected:
    //! An overridable function to set the colour of element index hi
    virtual std::array<float, 3> setColour (unsigned int hi)
    {
        std::array<float, 3> clr = { 0.0f, 0.0f, 0.0f };
        if (this->cm.getType() == morph::ColourMapType::Duochrome) {
            clr = this->cm.convert (this->dcolour[hi], this->dcolour2[hi]);
        } else {
            clr = this->cm.convert (this->dcolour[hi]);
        }
        return clr;
    }

    //! The retina data to visualize
    const retina<T>* ret;

    //! A copy of the scalarData which can be transformed suitably to be the z value of the surface
    std::vector<float> dcopy;
    //! A copy of the scalarData (or first field of vectorData), scaled to be a colour value
    std::vector<float> dcolour;
    std::vector<float> dcolour2;
};
