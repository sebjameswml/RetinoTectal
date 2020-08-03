/*
 * James et al model applied to Retinotectal projections.
 *
 * Author: Seb James <seb.james@sheffield.ac.uk>
 *
 * Date: Feb 2020.
 */

#ifndef FLT
// Check CMakeLists.txt to change to double or float
# error "Please define FLT when compiling (hint: See CMakeLists.txt)"
#endif

#include <iostream>
#define DBGSTREAM std::cout
#define DEBUG 1
#include <morph/MorphDbg.h>

//! A global variable used as a finished-the-loop flag. Global so that signalHandler can
//! set it true to break out of the loop early, but with the program completing correctly.
namespace sighandling {
    bool finished = false;
    bool user_interrupt = false;

    //! Signal handler to catch Ctrl-C
    void handler (int signum) {
        std::cerr << "User has interrupted the simulation. Finish up, save logs then exit...\n";
        //! Set true to finish early
        finished = true;
        user_interrupt = true;
    }
}

#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <csignal>
#include <chrono>

//! Our Retinotectal reaction diffusion class
#ifdef AXONCOMP
# include "rd_rettec.h"
#else
# include "rd_rettec_nocomp.h"
#endif

#include "morph/Vector.h"
// Shape analysis utilities
#include "morph/ShapeAnalysis.h"
//! Included for directory manipulation code
#include <morph/tools.h>
//! A jsoncpp-wrapping class for configuration.
#include <morph/Config.h>

#ifdef COMPILE_PLOTTING
# include <GLFW/glfw3.h>
// Scaling of values to be suitable for plotting
# include "morph/Scale.h"
// Colour maps!
# include "morph/ColourMap.h"
// A Visual gives you a 'visual scene'
# include "morph/Visual.h"
// All the Visual models here derive from VisualDataModel
# include "morph/VisualDataModel.h"
// Alias VisualDataModel<FLT>* as VdmPtr, to neaten code
typedef morph::VisualDataModel<FLT>* VdmPtr;
// We're visualizing HexGrids...
# include "morph/HexGridVisual.h"
// and doing quiver plots...
# include "morph/QuiverVisual.h"
// and a scatter plot.
# include "morph/ScatterVisual.h"

# include <morph/MathAlgo.h>

//! Helper function to save PNG images
void savePngs (const std::string& logpath, const std::string& name,
               unsigned int frameN, morph::Visual& v) {
    std::stringstream ff1;
    ff1 << logpath << "/" << name << "_";
    ff1 << std::setw(5) << std::setfill('0') << frameN;
    ff1 << ".png";
    v.saveImage (ff1.str());
}

//! Take the first element of the array and create a vector<vector<FLT>> to plot
std::vector<std::vector<FLT> > separateVectorField (std::vector<std::array<std::vector<FLT>, 2> >& f,
                                          unsigned int arrayIdx) {
    std::vector<std::vector<FLT> > vf;
    for (std::array<std::vector<FLT>, 2> fia : f) {
        std::vector<FLT> tmpv = fia[arrayIdx];
        vf.push_back (tmpv);
    }
    return vf;
}
#endif // COMPILE_PLOTTING

int main (int argc, char **argv)
{
    // register signal handler
    signal (SIGINT, sighandling::handler);
    signal (SIGTERM, sighandling::handler);

    // Randomly set the RNG seed
    srand (morph::Tools::randomSeed());

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " /path/to/params.json [/path/to/logdir]" << std::endl;
        return 1;
    }
    std::string paramsfile (argv[1]);

    // Set up a morph::Config object for reading configuration
    morph::Config conf(paramsfile);
    if (!conf.ready) {
        std::cerr << "Error setting up JSON config: " << conf.emsg << std::endl;
        return 1;
    }

    /*
     * Get simulation-wide parameters from JSON
     */
    const unsigned int steps = conf.getUInt ("steps", 1000UL);
    if (steps == 0) {
        std::cerr << "Not much point simulating 0 steps! Exiting." << std::endl;
        return 1;
    }
    // If logevery is 0, then log nothing to HDF5
    const unsigned int logevery = conf.getUInt ("logevery", 100UL);
    // Only start logging after simulation has got to this step:
    const unsigned int logfrom = conf.getUInt ("logfrom", 0UL);
    const float hextohex_d = conf.getFloat ("hextohex_d", 0.01f);
    const float boundaryFalloffDist = conf.getFloat ("boundaryFalloffDist", 0.01f);
    const std::string svgpath = conf.getString ("svgpath", "");
    // If svgpath is empty, then compute an elliptical boundary:
    const float ellipse_a = conf.getFloat ("ellipse_a", 1.0f);
    const float ellipse_b = conf.getFloat ("ellipse_b", 0.7f);
    // Full origin retina, or just a pie slice?
    const double ret_startangle = conf.getDouble ("ret_startangle", 0.0);
    const double ret_endangle = conf.getDouble ("ret_endangle", morph::TWO_PI_D);
    bool overwrite_logs = conf.getBool ("overwrite_logs", false);
    std::string logpath = conf.getString ("logpath", "fromfilename");
    std::string logbase = "";
    if (logpath == "fromfilename") {
        // Using json filename as logpath
        std::string justfile = paramsfile;
        // Remove trailing .json and leading directories
        std::vector<std::string> pth = morph::Tools::stringToVector (justfile, "/");
        justfile = pth.back();
        morph::Tools::searchReplace (".json", "", justfile);
        // Use logbase as the subdirectory into which this should go
        logbase = conf.getString ("logbase", "logs/");
        if (logbase.back() != '/') {
            logbase += '/';
        }
        logpath = logbase + justfile;
    }
    if (argc == 3) {
        std::string argpath(argv[2]);
        std::cerr << "Overriding the config-given logpath " << logpath << " with " << argpath << std::endl;
        logpath = argpath;
        if (overwrite_logs == true) {
            std::cerr << "WARNING: You set a command line log path.\n"
                      << "       : Note that the parameters config permits the program to OVERWRITE LOG\n"
                      << "       : FILES on each run (\"overwrite_logs\" is set to true)." << std::endl;
        }
    }

    // Used to initialise a
    const double aNoiseGain = conf.getDouble ("aNoiseGain", 0.1);
    const double aInitialOffset = conf.getDouble ("aInitialOffset", 0.1);
    const FLT dt = static_cast<FLT>(conf.getDouble ("dt", 0.00001));
    const FLT contour_threshold = conf.getDouble ("contour_threshold", 0.6);
    const double D = conf.getDouble ("D", 0.1);
    const double G = conf.getDouble ("G", 1.0);
    // sigma_gamma is the noise in determining gamma params from retinal neuron
    // positions. Applied once only at simulation setup.
    const double sigma_gamma = conf.getDouble ("sigma_gamma", 0.0);
    // sigma_rho is the noise introduced into the determination of the signalling
    // gradient. Applied on every simulation step.
    const double sigma_rho = conf.getDouble ("sigma_rho", 0.0);
    const FLT k = conf.getDouble ("k", 3.0);
    //const FLT l = conf.getDouble ("l", 1.0);
    //const FLT m = conf.getDouble ("m", 1e-8);

    DBG ("steps to simulate: " << steps);

    // Guidance molecule array of parameters:
    const Json::Value guid = conf.getArray("guidance");
    unsigned int M_GUID = static_cast<unsigned int>(guid.size());

#ifdef COMPILE_PLOTTING
    // Parameters from the config that apply only to plotting:
    const unsigned int plotevery = conf.getUInt ("plotevery", 10UL);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", false);
    unsigned int framecount = 0;
    // Which windows to plot?
    const bool plot_guide = conf.getBool ("plot_guide", true);
    const bool plot_contours = conf.getBool ("plot_contours", true);
    const bool plot_a_contours = conf.getBool ("plot_a_contours", true);
    const bool plot_a = conf.getBool ("plot_a", true);
    const bool plot_c = conf.getBool ("plot_c", true);
# ifndef AXONCOMP
    const bool plot_f = conf.getBool ("plot_f", false);
# endif
    const bool plot_n = conf.getBool ("plot_n", true);
    const bool plot_dr = conf.getBool ("plot_dr", true);
    const bool plot_guidegrad = conf.getBool ("plot_guidegrad", false);

    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height_default = static_cast<unsigned int>(0.8824f * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);

    // Set up the morph::Visual object
    morph::Visual v1 (win_width, win_height, "Retino-tectal simulation");
    v1.zNear = 0.001;
    v1.zFar = 50;
    v1.fov = 45;
    // If this is true, then mouse movements won't move the scene around
    v1.sceneLocked = conf.getBool ("sceneLocked", false);
    // Whole-scene offsetting
    v1.setZDefault (conf.getFloat ("z_default", -5.0f));
    v1.setSceneTransXY (conf.getFloat ("x_default", 0.0f), conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v1.scenetrans_stepsize = 0.5;
#endif // COMPILE_PLOTTING

    // Instantiate and set up the model object
#ifdef AXONCOMP
    RD_RetTec<FLT> RD;
#else
    RD_RetTec_NoComp<FLT> RD;
#endif
    RD.svgpath = svgpath;
    RD.ellipse_a = ellipse_a;
    RD.ellipse_b = ellipse_b;

    RD.ret_startangle = static_cast<FLT>(ret_startangle);
    RD.ret_endangle = static_cast<FLT>(ret_endangle);

    RD.logpath = logpath;
    // NB: Set .N, .M BEFORE RD.allocate().
    RD.N = conf.getUInt ("N", 0); // Number of RT axons is N
    RD.M = M_GUID; // Number of guidance molecules that are sculpted
    // Set up timestep
    RD.set_dt (dt);
    // Control the size of the hexes, and therefore the number of hexes in the grid
    RD.hextohex_d = hextohex_d;
    //RD.hexspan = hexspan; // auto-set in rd_barrel.h
    // Boundary fall-off distance
    RD.boundaryFalloffDist = boundaryFalloffDist;
    RD.aNoiseGain = aNoiseGain;
    RD.aInitialOffset = aInitialOffset;
    // Guidance molecule noise
    RD.mNoiseGain = conf.getDouble ("mNoiseGain", 0.0);
    RD.mNoiseSigma = conf.getDouble ("mNoiseSigma", 0.09);
    // After setting N and M, we can set up all the vectors in RD:
    RD.allocate();
    // After allocate(), we can set up the simulation parameters:
    RD.set_D (D);
    RD.G = G; // "gamma gain"
    RD.sigma_gamma = sigma_gamma;
    RD.sigma_rho = sigma_rho;
    RD.contour_threshold = contour_threshold;
    RD.k = k;
    RD.alpha_ = conf.getDouble ("alpha", 3.0);
    RD.beta_ = conf.getDouble ("beta", 20.0);
#ifdef AXONCOMP
    RD.epsilon = conf.getDouble ("epsilon", 0.2);
#else
    RD.s = conf.getDouble ("s", 1.0);
    RD.set_w (conf.getDouble ("w", 1.0));
    DBG ("Set RD.s to " << RD.s << " and RD.w to " << RD.get_w());
#endif

    // Index through guidance molecule parameters:
    for (unsigned int j = 0; j < guid.size(); ++j) {
        Json::Value v = guid[j];
        // What guidance molecule method will we use?
        std::string rmeth = v.get ("shape", "Sigmoid1D").asString();
        DBG ("guidance molecule shape: " << rmeth);
        if (rmeth == "Sigmoid1D") {
            RD.rhoMethod[j] = FieldShape::Sigmoid1D;
        } else if (rmeth == "Linear1D") {
            RD.rhoMethod[j] = FieldShape::Linear1D;
        } else if (rmeth == "Exponential1D") {
            RD.rhoMethod[j] = FieldShape::Exponential1D;
        } else if (rmeth == "Gauss1D") {
            RD.rhoMethod[j] = FieldShape::Gauss1D;
        } else if (rmeth == "Gauss2D") {
            RD.rhoMethod[j] = FieldShape::Gauss2D;
        } else if (rmeth == "CircLinear2D") {
            RD.rhoMethod[j] = FieldShape::CircLinear2D;
        }
        // Set up guidance molecule method parameters
        RD.guidance_gain.push_back (v.get("gain", 1.0).asDouble());
        DBG2 ("guidance modelecule gain: " << RD.guidance_gain.back());
        RD.guidance_phi.push_back (v.get("phi", 1.0).asDouble());
        RD.guidance_width.push_back (v.get("width", 1.0).asDouble());
        RD.guidance_offset.push_back (v.get("offset", 1.0).asDouble());
        RD.guidance_time_onset.push_back (v.get("time_onset", 0).asUInt());
    }

    // Now have the guidance molecule densities/gradients computed, call init()
    RD.init();

    // Create a log/png directory if necessary, and exit on any failures.
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            std::cerr << "Failed to create the logpath directory "
                      << logpath << " which does not exist." << std::endl;
            return 1;
        }
    } else {
        // Directory DOES exist. See if it contains a previous run and
        // exit without overwriting to avoid confusion.
        if (overwrite_logs == false
            && (morph::Tools::fileExists (logpath + "/params.json") == true
                || morph::Tools::fileExists (logpath + "/guidance.h5") == true
                || morph::Tools::fileExists (logpath + "/positions.h5") == true)) {
            std::cerr << "Seems like a previous simulation was logged in " << logpath << ".\n"
                      << "Please clean it out manually, choose another directory or set\n"
                      << "overwrite_logs to true in your parameters config JSON file." << std::endl;
            return 1;
        }
    }

    if (logevery > 0) {
        // As RD.allocate() as been called (and log directory has been
        // created/verified ready), positions, grid and guidance can be saved to file.
        RD.savePositions();
        RD.saveHG();
        RD.saveGuidance();
    }

#ifdef COMPILE_PLOTTING
    // Data scaling parameters
    float _m = 0.2;
    float _c = 0.0;
    // Z position scaling - how hilly/bumpy the visual will be.
    morph::Scale<FLT> zscale; zscale.setParams (_m/10.0f, _c/10.0f);
    // The second is the colour scaling.
    morph::Scale<FLT> cscale; cscale.setParams (_m, _c);
# ifndef AXONCOMP
    // scaling2
    morph::Scale<FLT> zscale2; zscale2.setParams (1.0f/5.0f, 0.0f);
    morph::Scale<FLT> cscale2; cscale2.setParams (1.0f, 0.0f);
# endif

    // Identifiers for the various VisualModels that will be added to the Visual scene
    unsigned int c_ctr_grid = 0;
    unsigned int a_ctr_grid = 0;
    unsigned int dr_grid = 0;
    unsigned int quiv_grid = 0;

    std::vector<unsigned int> guide_grids;
    std::vector<unsigned int> guidegrad_grids;

    // Spatial offset
    morph::Vector<FLT, 3> spatOff;

    // Start at a negative value which is determined by plot_a, plot_c and N.
    float xzero = 0.0f;
    xzero -= (plot_a == true ?  RD.hg->width() * sqrt(RD.N) : 0.0f);
    xzero -= (plot_c == true ?  RD.hg->width() * sqrt(RD.N) : 0.0f);
# ifndef AXONCOMP
    xzero -= (plot_f == true ?  RD.hg->width() * sqrt(RD.N) : 0.0f);
# endif

    // The a variable
    std::vector<unsigned int> agrids;
    unsigned int side = static_cast<unsigned int>(floor (sqrt (RD.N)));
    if (plot_a) {
        spatOff = {xzero, 0.0, 0.0 };
        for (unsigned int i = 0; i<RD.N; ++i) {
            spatOff[0] = xzero + RD.hg->width() * (i/side);
            spatOff[1] = RD.hg->width() * (i%side);
            unsigned int idx = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                                 RD.hg,
                                                                                 spatOff,
                                                                                 &RD.a[i],
                                                                                 zscale,
                                                                                 cscale,
                                                                                 morph::ColourMapType::Monochrome,
                                                                                 (float)i/(float)RD.N)); // hue
            agrids.push_back (idx);
        }
        xzero = spatOff[0] + RD.hg->width();
    }

    // The c variable
    std::vector<unsigned int> cgrids;
    if (plot_c) {
        spatOff = {xzero, 0.0, 0.0 };
        for (unsigned int i = 0; i<RD.N; ++i) {
            spatOff[0] = xzero + RD.hg->width() * (i/side);
            spatOff[1] = RD.hg->width() * (i%side);
            unsigned int idx = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                                 RD.hg,
                                                                                 spatOff,
                                                                                 &RD.c[i],
                                                                                 zscale,
                                                                                 cscale,
                                                                                 morph::ColourMapType::Monochrome,
                                                                                 (float)i/(float)RD.N)); // hue
            cgrids.push_back (idx);
        }
        xzero = spatOff[0] + RD.hg->width();
    }

# ifndef AXONCOMP
    // The f variable
    std::vector<unsigned int> fgrids;
    if (plot_f) {
        spatOff = {xzero, 0.0, 0.0 };
        for (unsigned int i = 0; i<RD.N; ++i) {
            spatOff[0] = xzero + RD.hg->width() * (i/side);
            spatOff[1] = RD.hg->width() * (i%side);
            unsigned int idx = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                                 RD.hg,
                                                                                 spatOff,
                                                                                 &RD.f[i],
                                                                                 zscale2,
                                                                                 cscale2,
                                                                                 morph::ColourMapType::Monochrome,
                                                                                 (float)i/(float)RD.N)); // hue
            fgrids.push_back (idx);
        }
        xzero = spatOff[0] + RD.hg->width();
    }
# endif
    // n
    unsigned int ngrid = 0;
    if (plot_n) {
        spatOff = { xzero, 0.0, 0.0 };
        ngrid = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                  RD.hg,
                                                                  spatOff,
                                                                  &RD.n,
                                                                  zscale,
                                                                  cscale,
                                                                  morph::ColourMapType::Plasma));
        xzero += RD.hg->width();
    }

    // Contours
    morph::Scale<FLT> null_zscale; null_zscale.setParams (0.0f, 0.0f);
    morph::Scale<FLT> ctr_cscale; ctr_cscale.setParams (1.0f, 0.0f);

    std::vector<FLT> zeromap (RD.nhex, static_cast<FLT>(0.0));

    std::vector<morph::Vector<FLT,3>> zerovecs;
    zerovecs.resize (RD.N);

    if (plot_contours) {
        spatOff = { xzero, 0.0, 0.0 };
        // special scaling for contours. flat in Z, but still colourful.
        // BUT, what I want is colours set by hue and i/N. That means a 'rainbow' colour map!
        c_ctr_grid = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                       RD.hg,
                                                                       spatOff,
                                                                       &zeromap,
                                                                       null_zscale,
                                                                       ctr_cscale,
                                                                       morph::ColourMapType::RainbowZeroBlack));
        xzero += RD.hg->width();
    }

    if (plot_a_contours) {
        spatOff = { xzero, 0.0, 0.0 };
        a_ctr_grid = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                       RD.hg,
                                                                       spatOff,
                                                                       &zeromap,
                                                                       null_zscale,
                                                                       ctr_cscale,
                                                                       morph::ColourMapType::RainbowZeroWhite));
        xzero += (1.2 * RD.hg->width());
    }

    if (plot_dr == true) {
        spatOff = { xzero, 0.0, 0.0 };
        dr_grid = v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                                    RD.hg,
                                                                    spatOff,
                                                                    &zeromap,
                                                                    null_zscale,
                                                                    ctr_cscale,
                                                                    morph::ColourMapType::RainbowZeroWhite));

        quiv_grid = v1.addVisualModel (new morph::QuiverVisual<FLT> (v1.shaderprog,
                                                                     &zerovecs,
                                                                     spatOff,
                                                                     &zerovecs,
                                                                     morph::ColourMapType::Fixed,
                                                                     0.2f));
        xzero +=  (1.2 * RD.hg->width());
    }

    // guidance expression
    if (plot_guide) {
        spatOff = { xzero, 0.0, 0.0 };
        morph::Scale<FLT> gd_cscale;
        gd_cscale.do_autoscale = true;
        // Plot gradients of the guidance effect g.
        for (unsigned int j = 0; j<RD.M; ++j) {
            v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                              RD.hg,
                                                              spatOff,
                                                              &RD.rho[j],
                                                              null_zscale,
                                                              gd_cscale,
                                                              morph::ColourMapType::Jet));
            spatOff[1] += 1.2f * RD.hg->depth();
        }

        // Plot coordinates of the Retinal neurons.
        xzero +=  (1.7 * RD.hg->width());
        spatOff = { xzero, 0.0, 0.0 };
        std::vector<morph::Vector<float, 3>> ret_coordinates;
        for (unsigned int c = 0; c < RD.ret_coords.size(); ++c) {
            std::array<float, 2> rc = RD.ret_coords[c];
            morph::Vector<float, 3> rc3;
            rc3[0] = rc[0];
            rc3[1] = rc[1];
            rc3[2] = 0.0f;
            ret_coordinates.push_back (rc3);
        }
        std::vector<float> neuronColourData;
        for (unsigned int i = 1; i <= RD.N; ++i) {
            neuronColourData.push_back ((float)i/(float)(RD.N+1));
        }
        float scatRad = RD.ring_d/10.0f;
        v1.addVisualModel (new morph::ScatterVisual<FLT> (v1.shaderprog, &ret_coordinates, spatOff,
                                                          &neuronColourData, scatRad, ctr_cscale,
                                                          morph::ColourMapType::RainbowZeroBlack));


        xzero +=  (1.2 * RD.hg->width());
    }

    // Now plot fields and redraw display
    if (plot_guidegrad) {
        spatOff = { xzero, 0.0, 0.0 };
        for (unsigned int j = 0; j<RD.M; ++j) {

            // gradient of guidance expression
            std::vector<std::vector<FLT> > gx = separateVectorField (RD.g[j], 0);
            std::vector<std::vector<FLT> > gy = separateVectorField (RD.g[j], 1);
            FLT ming = 1e7;
            FLT maxg = -1e7;
            if (plot_guidegrad) {
                // Determine scale of gx and gy so that a common scale can be
                // applied to both gradient_x and gradient_y.
                for (unsigned int hi=0; hi<RD.nhex; ++hi) {
                    morph::Hex* h = RD.hg->vhexen[hi];
                    if (h->onBoundary() == false) {
                        for (unsigned int i = 0; i<RD.N; ++i) {
                            if (gx[i][h->vi]>maxg) { maxg = gx[i][h->vi]; }
                            if (gx[i][h->vi]<ming) { ming = gx[i][h->vi]; }
                            if (gy[i][h->vi]>maxg) { maxg = gy[i][h->vi]; }
                            if (gy[i][h->vi]<ming) { ming = gy[i][h->vi]; }
                        }
                    }
                }
                DBG2 ("min g = " << ming << " and max g = " << maxg);
            }

            // Convert to a scaling object
            float gg_m, gg_c;
            gg_m = 1.0f/(float)(maxg-ming);
            gg_c = -(gg_m * ming);
            morph::Scale<float> gd_cscale; gd_cscale.setParams (gg_m, gg_c);
            // Create the grids
            v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                              RD.hg,
                                                              spatOff,
                                                              &gx[j],
                                                              null_zscale,
                                                              gd_cscale,
                                                              morph::ColourMapType::Jet));
            spatOff[0] += RD.hg->width();
            v1.addVisualModel (new morph::HexGridVisual<FLT> (v1.shaderprog,
                                                              RD.hg,
                                                              spatOff,
                                                              &gy[j],
                                                              null_zscale,
                                                              gd_cscale,
                                                              morph::ColourMapType::Jet));
            spatOff[0] -= RD.hg->width();
            spatOff[1] += RD.hg->depth();
        }
        xzero +=  (1.5 * RD.hg->width());
    }

    // Saving of t=0 images in log folder
    if ((RD.M > 0 && plot_guide) || plot_a) {
        savePngs (logpath, "sim", 0, v1);
    }

    // if using plotting, then set up the render clock
    std::chrono::steady_clock::time_point lastrender = std::chrono::steady_clock::now();

    // A pointer to access the data layer of HexGridVisual objects.
    morph::VisualDataModel<FLT>* mdlptr = (VdmPtr)0;
#endif // COMPILE_PLOTTING

    // Innocent until proven guilty
    sighandling::user_interrupt = false;
    conf.set ("crashed", false);

    try {
        // Start the loop
        sighandling::finished = false;
        while (sighandling::finished == false) {
            // Step the model
            RD.step();
            if ((RD.stepCount % 1000) == 0) {
                std::cout << RD.stepCount << " steps...\n";
            }

#ifdef COMPILE_PLOTTING
            if ((RD.stepCount % plotevery) == 0) {
                DBG2("Plot at step " << RD.stepCount);
                // Do a plot of the ctrs as found.
                std::vector<FLT> ctrmap = morph::ShapeAnalysis<FLT>::get_contour_map_nozero (RD.hg, RD.c, RD.contour_threshold);

                if (plot_contours) {
                    mdlptr = (VdmPtr)v1.getVisualModel (c_ctr_grid);
                    mdlptr->updateData (&ctrmap);
                }

                if (plot_a_contours) {
                    std::vector<FLT> actrmap = morph::ShapeAnalysis<FLT>::get_contour_map_nozero (RD.hg, RD.a, RD.contour_threshold);
                    mdlptr = (VdmPtr)v1.getVisualModel (a_ctr_grid);
                    mdlptr->updateData (&actrmap);
                }

                if (plot_a) {
                    for (unsigned int i = 0; i<RD.N; ++i) {
                        mdlptr = (VdmPtr)v1.getVisualModel (agrids[i]);
                        mdlptr->updateData (&RD.a[i]);
                    }
                }
                if (plot_c) {
                    for (unsigned int i = 0; i<RD.N; ++i) {
                        mdlptr = (VdmPtr)v1.getVisualModel (cgrids[i]);
                        mdlptr->updateData (&RD.c[i]);
                    }
                }
#ifndef AXONCOMP
                if (plot_f) {
                    for (unsigned int i = 0; i<RD.N; ++i) {
                        mdlptr = (VdmPtr)v1.getVisualModel (fgrids[i]);
                        mdlptr->updateData (&RD.f[i]);
                    }
                }
#endif
                if (plot_n) {
                    mdlptr = (VdmPtr)v1.getVisualModel (ngrid);
                    mdlptr->updateData (&RD.n);
                }
                if (plot_dr) {
                    RD.spatialAnalysis();
                    mdlptr = (VdmPtr)v1.getVisualModel (dr_grid);
                    mdlptr->updateData (&RD.regions);
                    // Plot the difference vectors here.
                    std::vector<morph::Vector<float, 3>> regcs;
                    for (auto rc : RD.reg_centroids) {
                        regcs.push_back ({rc.second.first, rc.second.second, 0.0f});
                    }
                    mdlptr = (VdmPtr)v1.getVisualModel (quiv_grid);
                    mdlptr->updateData (&regcs, &RD.tec_offsets);

                }
                // Save to PNG
                if (vidframes) {
                    savePngs (logpath, "sim", framecount, v1);
                    ++framecount;
                } else {
                    savePngs (logpath, "sim", RD.stepCount, v1);
                }
            }

            // rendering the graphics.
            std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
            if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
                glfwPollEvents();
                v1.render();
                lastrender = std::chrono::steady_clock::now();
            }
#endif // COMPILE_PLOTTING

            // Save data every 'logevery' steps
            if (logevery != 0 && RD.stepCount >= logfrom && (RD.stepCount == 1 || (RD.stepCount % logevery) == 0)) {
                DBG ("Logging data at step " << RD.stepCount);
                RD.save();
                // If spatial analysis, then add line here to do it
                RD.spatialAnalysis();
                // And save it
                RD.saveSpatial();
            }

            if (RD.stepCount > steps) {
                sighandling::finished = true;
            }
        }
    } catch (const std::exception& e) {
        // Set some stuff in the config as to what happened, so it'll get saved into params.conf
        conf.set ("crashed", true);
        std::stringstream ee;
        ee << e.what();
        conf.set ("exception_message", ee.str());
    }

    // Save out the sums.
    if (logevery > 0) {
        RD.savesums();
    }

    // Before saving the json, we'll place any additional useful info
    // in there, such as the FLT. If float_width is 4, then
    // results were computed with single precision, if 8, then double
    // precision was used. Also save various parameters from the RD system.
    conf.set ("float_width", (unsigned int)sizeof(FLT));
    std::string tnow = morph::Tools::timeNow();
    conf.set ("sim_ran_at_time", tnow.substr(0,tnow.size()-1));
    conf.set ("final_step", RD.stepCount);
    conf.set ("user_interrupt", sighandling::user_interrupt);
    conf.set ("hextohex_d", RD.hextohex_d);
    conf.set ("D", RD.get_D());
    conf.set ("k", RD.k);
    conf.set ("dt", RD.get_dt());
    // Call our function to place git information into root.
    //conf.insertGitInfo ("sim/");
    // Store the binary name and command argument into root, too.
    if (argc > 0) { conf.set("argv0", argv[0]); }
    if (argc > 1) { conf.set("argv1", argv[1]); }

    // We'll save a copy of the parameters for the simulation in the log directory as params.json
    const std::string paramsCopy = logpath + "/params.json";
    conf.write (paramsCopy);
    if (conf.ready == false) {
        std::cerr << "Warning: Something went wrong writing a copy of the params.json: " << conf.emsg << std::endl;
    }

#if 0
    // Extract contours
    std::vector<std::list<morph::Hex> > ctrs = morph::ShapeAnalysis<FLT>::get_contours (RD.hg, RD.c, RD.contour_threshold);
    {
        // Write each contour to a contours.h5 file
        std::stringstream ctrname;
        ctrname << logpath << "/contours.h5";
        morph::HdfData ctrdata(ctrname.str());
        unsigned int nctrs = ctrs.size();
        ctrdata.add_val ("/num_contours", nctrs);
        for (unsigned int ci = 0; ci < nctrs; ++ci) {
            std::vector<FLT> vx, vy;
            auto hi = ctrs[ci].begin();
            while (hi != ctrs[ci].end()) {
                vx.push_back (hi->x);
                vy.push_back (hi->y);
                ++hi;
            }
            std::stringstream ciss;
            ciss << ci;
            std::string pth = "/x" + ciss.str();
            ctrdata.add_contained_vals (pth.c_str(), vx);
            pth[1] = 'y';
            ctrdata.add_contained_vals (pth.c_str(), vy);

            // Generate hex grids from contours to obtain the size of the region enclosed by the contour
            morph::HexGrid* hg1 = new morph::HexGrid (RD.hextohex_d, RD.hexspan, 0, morph::HexDomainShape::Boundary);
            hg1->setBoundary (ctrs[ci]);
            pth[1] = 'n';
            ctrdata.add_val(pth.c_str(), hg1->num());
            delete hg1;
        }

        // Also extract the boundary of the main, enclosing hexgrid and write that.
        std::list<morph::Hex> outerBoundary = RD.hg->getBoundary();
        std::vector<FLT> vx, vy;
        auto bi = outerBoundary.begin();
        while (bi != outerBoundary.end()) {
            vx.push_back (bi->x);
            vy.push_back (bi->y);
            ++bi;
        }
        ctrdata.add_contained_vals ("/xb", vx);
        ctrdata.add_contained_vals ("/yb", vy);
    }
#endif

#ifdef COMPILE_PLOTTING
    if (sighandling::user_interrupt == false) {
        std::cout << "Press x in graphics window to exit.\n";
        v1.keepOpen();
    }
#endif

    return 0;
};
