
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"
#include "shapes/sphere.h"

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
	bool ThickLens(int start, int end, float &P1, float &P2, float &F1, float &F2) const;
  
private:
	// RealisticCamera Private Methods
	int nlens;
	
	// fundamental properties (assume at most 100 lens)
	float radius[100], axpos[100], lensN[100], aperture[100];
	float zcenter[100], zbound[100];
	
	float dist2film; // last lens to film, or called Z in the paper
	float exitpupil_r, exitpupil_z, exitpupil_C;
	
	Transform Raster2Camera;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H


