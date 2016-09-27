
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"

void TL_toImage(float P1, float P2, float F1, float F2, float &xt, float &yt, float &zt){
	float oldz = zt;
	float newz = P2 - 1.f / (1.f / (P2 - F2) - 1.f / (oldz - P1));
	float scale = (newz - P2) / (oldz - P1);

	xt *= scale, yt *= scale, zt = newz;
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
		float hither, float yon, 
		float sopen, float sclose, 
		float filmdistance, float aperture_diameter, string specfile, 
		float filmdiag, Film *f)
: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.
	FILE* camera_spec = fopen(specfile.c_str(), "r");

	char buf[300];
	while(fgets(buf, 300, camera_spec))
		if(buf[0] != '#') break;

	nlens = 0;
	lensN[nlens] = 1; // IoR for outside of the first lens (Void)
	nlens ++;
	
	sscanf(buf, "%f%f%f%f", &radius[nlens], &axpos[nlens], &lensN[nlens], &aperture[nlens]);
	nlens ++;

	int aper_id;

	while(fgets(buf, 300, camera_spec)){
		sscanf(buf, "%f%f%f%f", &radius[nlens], &axpos[nlens], &lensN[nlens], &aperture[nlens]);
		if(radius[nlens] == 0){
			aperture[nlens] = aperture_diameter;
			lensN[nlens] = 1.f;

			aper_id = nlens;
		}
		nlens ++;
	}
	nlens --;

	fprintf(stderr, "Number of Lens: %d\n", nlens);
	fclose(camera_spec);

	// Create a lot of spheres + an aperture to do intersection
	// Assume top point of the first lens is at origin, so film is in negative z

	float zdisp = 0;
	for(int i = 1; i <= nlens; i++){
		fprintf(stderr, "z-axis disp: %f\n", zdisp);
		
		zbound[i] = zdisp;
		zcenter[i] = zdisp - radius[i];
		zdisp -= axpos[i];
		
		fprintf(stderr, "lens %d: %f %f %f %f -> %f\n", i,
				radius[i], axpos[i], lensN[i], aperture[i], zcenter[i]);
	}

	dist2film = filmdistance;
	float zfilm = zdisp - dist2film;

	float rasterdiag = sqrt(1.*film->xResolution * film->xResolution + 1.*film->yResolution * film->yResolution);
	Raster2Camera = Translate(Vector(film->xResolution * filmdiag / rasterdiag / 2.f, 
									-film->yResolution * filmdiag / rasterdiag / 2.f, zfilm))
				  * Scale(-filmdiag / rasterdiag, filmdiag / rasterdiag, 1.f);
	
	printf("Aperture: %d\n", aper_id);
	
	float xP1, xP2, xF1, xF2;
	if(ThickLens(aper_id + 1, nlens, xP1, xP2, xF1, xF2) == 0)
		fprintf(stderr, "FAIL\n");

	exitpupil_r = aperture[aper_id] / 2.f;
	float tmp = 0;
	exitpupil_z = zbound[aper_id];
	
	printf("Aperture Stop: %f %f\n", exitpupil_r, exitpupil_z);
	
	TL_toImage(xP1, xP2, xF1, xF2, exitpupil_r, tmp, exitpupil_z);
	exitpupil_C = M_PI * exitpupil_r * exitpupil_r / ((exitpupil_z - zfilm) * (exitpupil_z - zfilm));

	printf("Exit Pupil: %f %f\n", exitpupil_r, exitpupil_z);
}

bool IntersectSphere(Ray Incident, float &tHit, float zcenter, float radius){
	Ray ray = Incident;
	ray.o.z -= zcenter;

    // Compute quadratic sphere coefficients
    float A = ray.d.x*ray.d.x + ray.d.y*ray.d.y + ray.d.z*ray.d.z;
    float B = 2 * (ray.d.x*ray.o.x + ray.d.y*ray.o.y + ray.d.z*ray.o.z);
    float C = ray.o.x*ray.o.x + ray.o.y*ray.o.y +
              ray.o.z*ray.o.z - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
	
    tHit = t0;
    if (t0 < ray.mint || (ray.o.z + t0 * ray.d.z) * radius < 0) {
        tHit = t1;
        if (tHit > ray.maxt) return false;
    }
	return true;
}

bool RealisticCamera::ThickLens(int start, int end, float &P1, float &P2, float &F1, float &F2) const {
	Point pt;
	Vector normal;

	float apermin = 100000.f;
	for(int i = start; i <= end; i++)
		apermin = min(apermin, aperture[i]);

	Point Pstart(apermin / 10.f, 0, zbound[start] + fabs(radius[start]));
	Vector Dstart(0, 0, -1);

	Ray Incident = Ray(Pstart, Dstart, 0);

	for(int i = start; i <= end; i++){
		float aper_r = aperture[i] / 2.f;
		
		if(radius[i] == 0){ // I am aperture
			float tHit = (zcenter[i] - Incident.o.z) / Incident.d.z;
			if(tHit < 0) return 0;

			pt = Incident(tHit);
			float r2 = pt.x * pt.x + pt.y * pt.y;
			if(r2 > aper_r * aper_r) return 0;

			Incident.o = pt;
		}
		else {
			// My own sphere intersection
			float tHit;
			if(!IntersectSphere(Incident, tHit, zcenter[i], radius[i])) return 0;

			pt = Incident(tHit);
			const float r2 = pt.x * pt.x + pt.y * pt.y;
			if(r2 > aper_r * aper_r) return 0;
			
			normal = (pt - Point(0, 0, zcenter[i])) / radius[i]; // The negative of below
		
			float eta = lensN[i - 1] / lensN[i]; // We are in material i - 1 going to i
			float c1 = -Dot(Incident.d, normal);
			float x1 = 1 - eta * eta * (1 - c1 * c1);
			if(x1 < 0) return 0;
			float c2 = sqrt(x1);

			// Incident.d should be a unit vector
			Incident = Ray(pt, (eta * Incident.d) + ((eta * c1 - c2) * normal), 0, INFINITY);
		}
	}
	
	F2 = Incident.o.z + Incident.d.z * ((0 - Incident.o.x) / Incident.d.x);
	P2 = Incident.o.z + Incident.d.z * ((Pstart.x - Incident.o.x) / Incident.d.x);

	// ========

	Pstart = Point(apermin / 10.f, 0, zbound[end] - fabs(radius[end]));
	Dstart = Vector(0, 0, 1);

	Incident = Ray(Pstart, Dstart, 0);

	for(int i = end; i >= start; i--){
		float aper_r = aperture[i] / 2.f;
		
		if(radius[i] == 0){ // I am aperture
			float tHit = (zcenter[i] - Incident.o.z) / Incident.d.z;
			if(tHit < 0) return 0;

			pt = Incident(tHit);
			float r2 = pt.x * pt.x + pt.y * pt.y;
			if(r2 > aper_r * aper_r) return 0;

			Incident.o = pt;
		}
		else {
			// My own sphere intersection
			float tHit;
			if(!IntersectSphere(Incident, tHit, zcenter[i], radius[i])) return 0;

			pt = Incident(tHit);
			const float r2 = pt.x * pt.x + pt.y  *pt.y;
			if(r2 > aper_r * aper_r) return 0;
			
			normal = -(pt - Point(0, 0, zcenter[i])) / radius[i];
		
			float eta = lensN[i] / lensN[i - 1];
			float c1 = -Dot(Incident.d, normal);
			float x1 = 1 - eta * eta * (1 - c1 * c1);
			if(x1 < 0) return 0;
			float c2 = sqrt(x1);

			// Incident.d should be a unit vector
			Incident = Ray(pt, (eta * Incident.d) + ((eta * c1 - c2) * normal), 0, INFINITY);
		}
	}

	F1 = Incident.o.z + Incident.d.z * ((0 - Incident.o.x) / Incident.d.x);
	P1 = Incident.o.z + Incident.d.z * ((Pstart.x - Incident.o.x) / Incident.d.x);

	fprintf(stderr, "F1:%f P1:%f P2:%f F2:%f\n", F1, P1, P2, F2);

	return 1;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// YOUR CODE HERE -- make that ray!

	// use sample->imageX and sample->imageY to get raster-space coordinates
	// of the sample point on the film.
	// use sample->lensU and sample->lensV to get a sample position on the lens

	Point Pras(sample.imageX, sample.imageY, 0);
	Point Pstart = Raster2Camera(Pras);

	float lensU, lensV;
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);

	Vector Initdir = Point(lensU * exitpupil_r, lensV * exitpupil_r, exitpupil_z) - Pstart;
	Initdir = Normalize(Initdir);

	Ray Incident(Pstart, Initdir, 0);

	Point pt;
	Vector normal;

	for(int i = nlens; i >= 1; i--) {
		float aper_r = aperture[i] / 2.f;
		
		if(radius[i] == 0){ // I am aperture
			float tHit = (zcenter[i] - Incident.o.z) / Incident.d.z;
			if(tHit < 0) return 0;

			pt = Incident(tHit);
			float r2 = pt.x * pt.x + pt.y * pt.y;
			if(r2 > aper_r * aper_r) return 0;

			Incident.o = pt;
		}
		else {
			// My own sphere intersection
			float tHit;
			if(!IntersectSphere(Incident, tHit, zcenter[i], radius[i])) return 0;

			pt = Incident(tHit);
			const float r2 = pt.x * pt.x + pt.y  *pt.y;
			if(r2 > aper_r * aper_r) return 0;
			
			normal = -(pt - Point(0, 0, zcenter[i])) / radius[i];
		
			float eta = lensN[i] / lensN[i-1];
			float c1 = -Dot(Incident.d, normal);
			float x1 = 1 - eta * eta * (1 - c1 * c1);
			if(x1 < 0) return 0;
			float c2 = sqrt(x1);

			// Incident.d should be a unit vector
			Incident = Ray(pt, (eta * Incident.d) + ((eta * c1 - c2) * normal), 0, INFINITY);
		}
	}

	*ray = CameraToWorld(Incident);
	return exitpupil_C * (Initdir.z * Initdir.z * Initdir.z * Initdir.z);
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
		const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
			shutterclose != -1 && filmdistance!= -1);

	if (specfile == "") {
		Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
			shutteropen, shutterclose, filmdistance, fstop, 
			specfile, filmdiag, film);
}
