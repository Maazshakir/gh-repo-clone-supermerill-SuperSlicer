#include "../ClipperUtils.hpp"
#include "../PolylineCollection.hpp"
#include "../Surface.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

#include "FillGyroid.hpp"

namespace Slic3r {

// 0 = auto (from 2 to ~6). note: if not in auto, the density is not properly calculated anymore.
static int const staticNbLines = 0;


static Polyline makeLineVert(double xPos, double yPos, double width, double height, double currentXBegin, double segmentSize, coord_t scaleFactor, double zCs, double zSn, bool flip, double decal=0){
	
	double maxSlope = abs(abs(zCs)-abs(zSn));
	Polyline polyline;
	polyline.points.push_back(Point(coord_t((std::max(std::min(currentXBegin, xPos+width),xPos) + decal) * scaleFactor), coord_t(yPos * scaleFactor)));
	for(double y=yPos;y<yPos+height+segmentSize;y+=segmentSize){
		if(y>yPos+height) y = yPos+height;
		double ySn = sin(y +(zCs<0?3.14:0) + 3.14);
		double yCs = cos(y +(zCs<0?3.14:0) + 3.14+(!flip?0:3.14));

		double a = ySn;
		double b = -zCs;
		double res = zSn*yCs;
		double r = sqrt(a*a + b*b);
		double x = asin(a/r) + asin(res/r) +3.14;
		x += currentXBegin;

		if(x < xPos){
			// needNewLine= true;
			x = xPos;
		}
		if(x > xPos+width){
			// needNewLine= true;
			x = xPos+width;
		}
		{
			double ydeviation = (flip?-1:1)*(zSn>0?-1:1)*decal*(1-maxSlope)*(res/r - a/r);
			polyline.points.push_back(Point(coord_t((x+decal) * scaleFactor), coord_t((y + ydeviation) * scaleFactor)));
		}
	}
	
	return polyline;
}

static Polyline makeLineHori(double xPos, double yPos, double width, double height, double currentYBegin, double segmentSize, coord_t scaleFactor, double zCs, double zSn, bool flip, double decal=0){
	double maxSlope = abs(abs(zCs)-abs(zSn));
	Polyline polyline;
	polyline.points.push_back(Point(coord_t(xPos * scaleFactor), coord_t((std::max(std::min(currentYBegin, yPos+height),yPos)+decal) * scaleFactor)));
	for(double x=xPos;x<xPos+width+segmentSize;x+=segmentSize){
		if(x>xPos+width) x = xPos+width;
		double xSn = sin(x +(zSn<0?3.14:0) +(flip?0:3.14));
		double xCs = cos(x +(zSn<0?3.14:0) );
		double a = xCs;
		double b = -zSn;
		double res = zCs*xSn;
		double r = sqrt(a*a + b*b);
		double y = asin(a/r) + asin(res/r) +3.14/2;
		y += currentYBegin;
		double xdeviation = (flip?-1:1)*(zCs>0?-1:1)*decal*(1-maxSlope)*(res/r - a/r);
		polyline.points.push_back(Point(coord_t((x + xdeviation) * scaleFactor), coord_t((std::max(std::min(y, yPos+height),yPos)+decal) * scaleFactor)));
	}
	
	return polyline;
}

static inline void correctOrderAndAdd(const int num, Polyline &poly, Polylines &array){
	if(num%2==0){
		Points temp(poly.points.rbegin(), poly.points.rend());
		poly.points.assign(temp.begin(),temp.end());
	}
	array.push_back(poly);
}

// Generate a set of curves (array of array of 2d points) that describe a
// horizontal slice of a truncated regular octahedron with a specified
// grid square size.
static Polylines makeGrid(coord_t gridZ, double density, double layer_width, double layer_height, size_t gridWidth, size_t gridHeight, size_t curveType)
{
	//reduce the density a bit to have a 90% infill not over-extruding.
	density = density * 0.8;
	
    coord_t  scaleFactor = coord_t(scale_(layer_width) / density);;
    Polylines result;
	Polyline *polyline2;
	double segmentSize = density/2;
	double decal = layer_width*density;
	double xPos = 0, yPos=0, width=gridWidth, height=gridHeight;
	 //scale factor for 5% : 8 712 388
	 // 1z = 10^-6 mm ?
	double z = gridZ/(1.0 * scaleFactor);
	// std::cout<<"gridZ= "<<gridZ<<", z? "<<z<<", scaleFactor= "<<scaleFactor<<", segmentSize= "<<segmentSize<<", decal= "<<decal<<" scale:"<<scale_(1)<<std::endl;
	double zSn = sin(z);
	double zCs = cos(z);
	double maxSlope = abs(abs(zCs)-abs(zSn));
	//be sure to don't be too near each other)
	while(maxSlope<0.05){
		// move away
		bool sameSign = zCs*zSn>=0;
		bool absCosGtAbsSin = abs(zCs)>abs(zSn);
		if(sameSign && !absCosGtAbsSin || !sameSign && absCosGtAbsSin){
			// then add
			z += 0.001;
		}else{
			// then remove
			z -= 0.001;
		}
		zSn = sin(z);
		zCs = cos(z);
		maxSlope = abs(abs(zCs)-abs(zSn));
		for(int i=0;i<100000000;i++){i%2==0?density+=0.001:density-=0.001;}
	}
	

	int numLine = 0;
	//nbLines depends of layer_width, layer_height, density and maxSlope 
	int nbLines = staticNbLines;
	if(staticNbLines==0){
		if(density>0.6){
			nbLines = 2;
		}else{
			//compute the coverage of the current layer
			double nbLineMore = 1 + (1 - maxSlope);
			nbLineMore *= nbLineMore;
			nbLineMore *= nbLineMore;
			nbLineMore = nbLineMore -1;
			nbLineMore *= 0.1 * layer_height / (layer_width * sqrt(density));
			nbLines = 2 + nbLineMore;
			// std::cout<<"|maxSlope="<<maxSlope<<"; nbLineMore= "<<nbLineMore<<std::endl;
		}
	}
	
	if(abs(zSn)<=abs(zCs)){
		//vertical
		//begin to first one
		int iter = 1;
		double currentXBegin = xPos - PI/2;
		currentXBegin = PI*(int)(currentXBegin/PI -1);
		iter = (int)(currentXBegin/PI +1)%2;
		bool flip = iter%2==1;
		// bool needNewLine =false;
		while(currentXBegin<xPos+width-PI/2){
			
			for(int num=0;num<nbLines;num++){
				double movex = -decal * (nbLines-1) + num*decal*2;
				correctOrderAndAdd(numLine, makeLineVert(xPos, yPos, width, height, currentXBegin, segmentSize, scaleFactor, zCs, zSn, flip, movex), result);
				numLine++;
			}
			
			//then, return by the other side
			iter++;
			currentXBegin = currentXBegin + PI;
			flip = iter%2==1;
			
			if(currentXBegin < xPos+width-PI/2){
				
				for(int num=0;num<nbLines;num++){
					double movex = -decal * (nbLines-1) + num*decal*2;
					correctOrderAndAdd(numLine, makeLineVert(xPos, yPos, width, height, currentXBegin, segmentSize, scaleFactor, zCs, zSn, flip, movex), result);
					numLine++;
				}

				// relance
				iter++;
				currentXBegin = currentXBegin + PI;
				flip = iter%2==1;
			}
		}
	}else{
		//horizontal
		

		//begin to first one
		int iter = 1;
		//search first line output
		double currentYBegin = yPos ;
		currentYBegin = PI*(int)(currentYBegin/PI -0);
		iter = (int)(currentYBegin/PI +1)%2;
		
		bool flip = iter%2==1;
		
		
		while(currentYBegin < yPos+width){

			for(int num=0;num<nbLines;num++){
				double movex = -decal * (nbLines-1) + num*decal*2;
				correctOrderAndAdd(numLine, makeLineHori(xPos, yPos, width, height, currentYBegin, segmentSize, scaleFactor, zCs, zSn, flip, movex), result);
				numLine++;
			}
		
			//then, return by the other side
			iter++;
			currentYBegin = currentYBegin + PI;
			flip = iter%2==1;
			
			if(currentYBegin<yPos+width){
				
				for(int num=0;num<nbLines;num++){
					double movex = -decal * (nbLines-1) + num*decal*2;
					correctOrderAndAdd(numLine, makeLineHori(xPos, yPos, width, height, currentYBegin, segmentSize, scaleFactor, zCs, zSn, flip, movex), result);
					numLine++;
				}
				
				//relance
				iter++;
				currentYBegin = currentYBegin + PI;
				flip = iter%2==1;
			}
		}
	}
	
    return result;
}

void FillGyroid::_fill_surface_single(
    const FillParams                &params, 
    unsigned int                     thickness_layers,
    const std::pair<float, Point>   &direction, 
    ExPolygon                       &expolygon, 
    Polylines                       &polylines_out)
{
    // no rotation is supported for this infill pattern
    BoundingBox bb = expolygon.contour.bounding_box();
    coord_t     distance = coord_t(scale_(this->spacing) / params.density);

	std::cout<<"thickness_layers="<<thickness_layers<<", spacing="<<this->spacing<<", scale_(this->spacing)="<<scale_(this->spacing)<<", params.density="<<params.density<<", distance="<<distance<<std::endl;
    // align bounding box to a multiple of our honeycomb grid module
    // (a module is 2*$distance since one $distance half-module is 
    // growing while the other $distance half-module is shrinking)
    bb.merge(_align_to_grid(bb.min, Point(2*distance, 2*distance)));
    
    // generate pattern
    Polylines   polylines = makeGrid(
        scale_(this->z),
        params.density,
        this->spacing,
		this->layer_height,
        ceil(bb.size().x / distance) + 1,
        ceil(bb.size().y / distance) + 1,
        ((this->layer_id/thickness_layers) % 2) + 1);
    
    // move pattern in place
    for (Polylines::iterator it = polylines.begin(); it != polylines.end(); ++ it)
        it->translate(bb.min.x, bb.min.y);

    // clip pattern to boundaries
    polylines = intersection_pl(polylines, (Polygons)expolygon);

    // connect lines
    if (! params.dont_connect && ! polylines.empty()) { // prevent calling leftmost_point() on empty collections
        ExPolygon expolygon_off;
        {
            ExPolygons expolygons_off = offset_ex(expolygon, SCALED_EPSILON);
            if (! expolygons_off.empty()) {
                // When expanding a polygon, the number of islands could only shrink. Therefore the offset_ex shall generate exactly one expanded island for one input island.
                assert(expolygons_off.size() == 1);
                std::swap(expolygon_off, expolygons_off.front());
            }
        }
        Polylines chained = PolylineCollection::chained_path_from(
#if SLIC3R_CPPVER >= 11
            std::move(polylines), 
#else
            polylines,
#endif
            PolylineCollection::leftmost_point(polylines), false); // reverse allowed
        bool first = true;
        for (Polylines::iterator it_polyline = chained.begin(); it_polyline != chained.end(); ++ it_polyline) {
            if (! first) {
                // Try to connect the lines.
                Points &pts_end = polylines_out.back().points;
                const Point &first_point = it_polyline->points.front();
                const Point &last_point = pts_end.back();
                // TODO: we should also check that both points are on a fill_boundary to avoid 
                // connecting paths on the boundaries of internal regions
				// question : why not using this->link_max_length?
                if (first_point.distance_to(last_point) <= 4 * distance && 
                    expolygon_off.contains(Line(last_point, first_point))) {
                    // Append the polyline.
                    pts_end.insert(pts_end.end(), it_polyline->points.begin(), it_polyline->points.end());
                    continue;
                }
            }
            // The lines cannot be connected.
#if SLIC3R_CPPVER >= 11
            polylines_out.push_back(std::move(*it_polyline));
#else
            polylines_out.push_back(Polyline());
            std::swap(polylines_out.back(), *it_polyline);
#endif
            first = false;
        }
    }
}

} // namespace Slic3r
