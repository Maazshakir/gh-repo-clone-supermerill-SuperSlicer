#ifndef slic3r_FillGyroid_hpp_
#define slic3r_FillGyroid_hpp_

#include <map>

#include "../libslic3r.h"

#include "FillBase.hpp"

namespace Slic3r {

class FillGyroid : public Fill
{
public:

	FillGyroid(){ wall_nb_lines = -1; scaling = 1.12f; }
    virtual Fill* clone() const { return new FillGyroid(*this); };
    virtual ~FillGyroid() {}

	// require bridge flow since most of this pattern hangs in air
    virtual bool use_bridge_flow() const { return true; }

protected:

	// 0 = auto (from 0 to ~4).
	//neg : auto whit min = -GYROID_WALL_NB_LINES
	//  note: if not -2, the density is not properly calculated anymore.
	int wall_nb_lines;
	
	float scaling;


	virtual void _fill_surface_single(
	    const FillParams                &params, 
	    unsigned int                     thickness_layers,
	    const std::pair<float, Point>   &direction, 
	    ExPolygon                       &expolygon, 
	    Polylines                       &polylines_out);
	
	// create the gyroid grid to clip.
	Polylines makeGrid(coord_t gridZ, double density, double layer_width, double layer_height, size_t gridWidth, size_t gridHeight, size_t curveType);
	//add line poly in reverse if needed into array
	inline void correctOrderAndAdd(const int num, Polyline &poly, Polylines &array);
	//create a curved horinzontal line  (for each x, compute y)
	Polyline makeLineHori(double xPos, double yPos, double width, double height, 
		double currentYBegin, double segmentSize, coord_t scaleFactor, 
		double zCs, double zSn, 
		bool flip, double decal=0);
	//create a curved vertival line (for each y, compute x)
	Polyline makeLineVert(double xPos, double yPos, double width, double height, 
		double currentXBegin, double segmentSize, coord_t scaleFactor, 
		double zCs, double zSn, 
		bool flip, double decal=0);

};


class FillGyroidThin : public FillGyroid
{
public:
	FillGyroidThin(){ wall_nb_lines = 1; scaling = 1.75; }
    virtual Fill* clone() const { return new FillGyroidThin(*this); };
    virtual ~FillGyroidThin() {}

};

class FillGyroidThick : public FillGyroid
{
public:
	FillGyroidThick(){ wall_nb_lines = -2; scaling = 0.59f; }
    virtual Fill* clone() const { return new FillGyroidThick(*this); };
    virtual ~FillGyroidThick() {}

};

} // namespace Slic3r

#endif // slic3r_FillGyroid_hpp_