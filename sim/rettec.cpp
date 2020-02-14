/*
 * James et al model applied to Retinotectal projections.
 *
 * Author: Seb James <seb.james@sheffield.ac.uk>
 *
 * Date: Feb 2020.
 */

/*!
 * This will be passed as the template argument for RD_Plot and RD and
 * should be defined when compiling.
 */
#ifndef FLT
// Check CMakeLists.txt to change to double or float
# error "Please define FLT when compiling (hint: See CMakeLists.txt)"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <chrono>
using namespace std::chrono;
using std::chrono::steady_clock;

//! Include the relevant reaction diffusion class
#include "rd_rettec.h"

#include "morph/ShapeAnalysis.h"
using morph::ShapeAnalysis;

#ifdef COMPILE_PLOTTING
# include <GLFW/glfw3.h>
# include "morph/Visual.h"
using morph::Visual;
# include <morph/MathAlgo.h>
using morph::MathAlgo;

//! Helper function to save PNG images
void savePngs (const std::string& logpath, const std::string& name,
               unsigned int frameN, Visual& v) {
    std::stringstream ff1;
    ff1 << logpath << "/" << name<< "_";
    ff1 << std::setw(5) << std::setfill('0') << frameN;
    ff1 << ".png";
    v.saveImage (ff1.str());
}

//! Take the first element of the array and create a vector<vector<FLT>> to plot
vector<vector<FLT> > separateVectorField (vector<array<vector<FLT>, 2> >& f,
                                          unsigned int arrayIdx) {
    vector<vector<FLT> > vf;
    for (array<vector<FLT>, 2> fia : f) {
        vector<FLT> tmpv = fia[arrayIdx];
        vf.push_back (tmpv);
    }
    return vf;
}
#endif

//! Included for directory manipulation code
#include "morph/tools.h"

//! A jsoncpp-wrapping class for configuration.
#include "morph/Config.h"

using namespace std;

int main (int argc, char **argv)
{
    // Randomly set the RNG seed
    srand (morph::Tools::randomSeed());

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " /path/to/params.json [/path/to/logdir]" << endl;
        return 1;
    }
    string paramsfile (argv[1]);

    // Set up a morph::Config object for reading configuration
    morph::Config conf(paramsfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
        return 1;
    }

    /*
     * Get simulation-wide parameters from JSON
     */
    const unsigned int steps = conf.getUInt ("steps", 1000UL);
    if (steps == 0) {
        cerr << "Not much point simulating 0 steps! Exiting." << endl;
        return 1;
    }
    const unsigned int logevery = conf.getUInt ("logevery", 100UL);
    if (logevery == 0) {
        cerr << "Can't log every 0 steps. Exiting." << endl;
        return 1;
    }
    const float hextohex_d = conf.getFloat ("hextohex_d", 0.01f);
    const float hexspan = conf.getFloat ("hexspan", 4.0f);
    const float boundaryFalloffDist = conf.getFloat ("boundaryFalloffDist", 0.01f);
    const string svgpath = conf.getString ("svgpath", "./ellipse.svg");
    bool overwrite_logs = conf.getBool ("overwrite_logs", false);
    string logpath = conf.getString ("logpath", "fromfilename");
    string logbase = "";
    if (logpath == "fromfilename") {
        // Using json filename as logpath
        string justfile = paramsfile;
        // Remove trailing .json and leading directories
        vector<string> pth = morph::Tools::stringToVector (justfile, "/");
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
        string argpath(argv[2]);
        cerr << "Overriding the config-given logpath " << logpath << " with " << argpath << endl;
        logpath = argpath;
        if (overwrite_logs == true) {
            cerr << "WARNING: You set a command line log path.\n"
                 << "       : Note that the parameters config permits the program to OVERWRITE LOG\n"
                 << "       : FILES on each run (\"overwrite_logs\" is set to true)." << endl;
        }
    }

    // Used to initialise a
    const double aNoiseGain = conf.getDouble ("aNoiseGain", 0.1);
    const double aInitialOffset = conf.getDouble ("aInitialOffset", 0.1);
    const FLT dt = static_cast<FLT>(conf.getDouble ("dt", 0.00001));
    const FLT contour_threshold = conf.getDouble ("contour_threshold", 0.6);
    const double D = conf.getDouble ("D", 0.1);
    const FLT k = conf.getDouble ("k", 3.0);
    const FLT l = conf.getDouble ("l", 1.0);
    const FLT m = conf.getDouble ("m", 1e-8);

    DBG ("steps to simulate: " << steps);

    // Retino-tectal projections. Same as number of hexes?
    const Json::Value rts = conf.getArray ("rt");
    unsigned int N_Axons =  conf.getUInt ("N", 0);
    if (N_Axons == 0) {
        cerr << "Zero retinotectal axons makes no sense for this simulation. Exiting."
             << endl;
        return 1;
    }

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
    const bool plot_n = conf.getBool ("plot_n", true);
    const bool plot_dr = conf.getBool ("plot_dr", true);
    const bool plot_guidegrad = conf.getBool ("plot_guidegrad", false);

    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height = static_cast<unsigned int>(0.8824f * (float)win_width);

    // Set up the morph::Visual object
    Visual plt (win_width, win_height, "Retino-tectal simulation");
    plt.zNear = 0.001;
    plt.zFar = 50;
    plt.fov = 45;
    plt.setZDefault (10.0);
#endif

    // Instantiate and set up the model object
    RD_RetTec<FLT> RD;
    RD.svgpath = svgpath;
    RD.logpath = logpath;
    // NB: Set .N, .M BEFORE RD.allocate().
    RD.N = N_Axons; // Number of RT axons.
    RD.M = M_GUID; // Number of guidance molecules that are sculpted
    // Set up timestep
    RD.set_dt (dt);
    // Control the size of the hexes, and therefore the number of hexes in the grid
    RD.hextohex_d = hextohex_d;
    RD.hexspan = hexspan;
    // Boundary fall-off distance
    RD.boundaryFalloffDist = boundaryFalloffDist;
    RD.aNoiseGain = aNoiseGain;
    RD.aInitialOffset = aInitialOffset;
    // After setting N and M, we can set up all the vectors in RD:
    RD.allocate();
    // After allocate(), we can set up parameters:
    RD.set_D (D);
    RD.l = l;
    RD.m = m;
    RD.E = static_cast<FLT>(0.0);
    RD.contour_threshold = contour_threshold;
    RD.k = k;

    // Index through thalamocortical fields, setting params:
    for (unsigned int i = 0; i < rts.size(); ++i) {
        Json::Value v = rts[i];
        RD.alpha[i] = v.get("alpha", 0.0).asDouble();
        RD.beta[i] = v.get("beta", 0.0).asDouble();

        // Sets up mask for initial branching density
        GaussParams<FLT> gp;
        gp.gain = v.get("gaininit", 1.0).asDouble();
        gp.sigma = v.get("sigmainit", 0.0).asDouble();
        gp.x = v.get("xinit", 0.0).asDouble();
        gp.y = v.get("yinit", 0.0).asDouble();
        RD.initmasks.push_back (gp);
        RD.epsilon[i] = v.get("epsilon", 0.0).asDouble();
        DBG2 ("Set RD.epsilon["<<i<<"] to " << RD.epsilon[i]);
    }

    // Index through guidance molecule parameters:
    for (unsigned int j = 0; j < guid.size(); ++j) {
        Json::Value v = guid[j];
        // What guidance molecule method will we use?
        string rmeth = v.get ("shape", "Sigmoid1D").asString();
        DBG2 ("guidance molecule shape: " << rmeth);
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

    // Which of the gammas is the "group" defining gamma?
    const unsigned int groupgamma = conf.getUInt ("groupgamma", 0UL);

    // Set up the interaction parameters between the different TC
    // populations and the guidance molecules (aka gamma).
    int paramRtn = 0;
    for (unsigned int i = 0; i < rts.size(); ++i) {
        Json::Value rtv = rts[i];
        Json::Value gamma = rtv["gamma"];
        Json::Value rtname = rtv["name"];
        for (unsigned int j = 0; j < guid.size(); ++j) {
            // Set up gamma values using a setter which checks we
            // don't set a value that's off the end of the gamma
            // container.
            DBG2 ("Set gamma for guidance " << j << " over axon " << i << " = " << gamma[j]);
            paramRtn += RD.setGamma (j, i, gamma[j].asDouble(), groupgamma);
        }
        // Make a map of name to float id value
        RD.rtnames[(FLT)i/(FLT)rts.size()] = rtname.asString();
    }

    if (paramRtn && M_GUID>0) {
        cerr << "Something went wrong setting gamma values" << endl;
        return paramRtn;
    }

    // Now have the guidance molecule densities and their gradients computed, call init()
    RD.init();

    /*
     * Now create a log directory if necessary, and exit on any
     * failures.
     */
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            cerr << "Failed to create the logpath directory "
                 << logpath << " which does not exist."<< endl;
            return 1;
        }
    } else {
        // Directory DOES exist. See if it contains a previous run and
        // exit without overwriting to avoid confusion.
        if (overwrite_logs == false
            && (morph::Tools::fileExists (logpath + "/params.json") == true
                || morph::Tools::fileExists (logpath + "/guidance.h5") == true
                || morph::Tools::fileExists (logpath + "/positions.h5") == true)) {
            cerr << "Seems like a previous simulation was logged in " << logpath << ".\n"
                 << "Please clean it out manually, choose another directory or set\n"
                 << "overwrite_logs to true in your parameters config JSON file." << endl;
            return 1;
        }
    }

    // As RD.allocate() as been called (and log directory has been
    // created/verified ready), positions can be saved to file.
    RD.savePositions();
    RD.saveHG();
    // Save the guidance molecules now.
    RD.saveGuidance();

#ifdef COMPILE_PLOTTING

    // Data scaling parameters
    float _m = 0.2;
    float _c = 0.0;
    const array<float, 4> scaling = { _m/10, _c/10, _m, _c };

    // HERE, add HexGridVisuals...
    unsigned int c_ctr_grid = 0; // one only
    unsigned int a_ctr_grid = 0;
#if 0
    unsigned int dr_grid = 0;
#endif
    vector<unsigned int> guide_grids;
    vector<unsigned int> guidegrad_grids;

    // Spatial offset
    array<float, 3> spatOff;
    float xzero = 0.0f;

    // The a variable
    vector<unsigned int> agrids;
    unsigned int side = static_cast<unsigned int>(floor (sqrt (RD.N)));
    if (plot_a) {
        spatOff = {xzero, 0.0, 0.0 };
        for (unsigned int i = 0; i<RD.N; ++i) {
            spatOff[0] = xzero + RD.hg->width() * (i/side);
            spatOff[1] = RD.hg->width() * (i%side);
            unsigned int idx = plt.addHexGridVisual (RD.hg, spatOff, RD.a[i], scaling);
            agrids.push_back (idx);
        }
        xzero = spatOff[0] + RD.hg->width();
    }

    // The c variable
    vector<unsigned int> cgrids;
    if (plot_c) {
        spatOff = {xzero, 0.0, 0.0 };
        for (unsigned int i = 0; i<RD.N; ++i) {
            spatOff[0] = xzero + RD.hg->width() * (i/side);
            spatOff[1] = RD.hg->width() * (i%side);
            unsigned int idx = plt.addHexGridVisual (RD.hg, spatOff, RD.c[i], scaling);
            cgrids.push_back (idx);
        }
        xzero = spatOff[0] + RD.hg->width();
    }

    // n
    unsigned int ngrid = 0;
    if (plot_n) {
        spatOff = { xzero, 0.0, 0.0 };
        ngrid = plt.addHexGridVisual (RD.hg, spatOff, RD.n, scaling);
        xzero += RD.hg->width();
    }

    // Contours
    const array<float, 4> ctr_scaling = { 0.0f, 0.0f, 1.0f, 0.0f };
    vector<FLT> zeromap (RD.nhex, static_cast<FLT>(0.0));
    if (plot_contours) {
        spatOff = { xzero, 0.0, 0.0 };
        // special scaling for contours. flat in Z, but still colourful
        c_ctr_grid = plt.addHexGridVisual (RD.hg, spatOff, zeromap, ctr_scaling);
        xzero += RD.hg->width();
    }

    if (plot_a_contours) {
        spatOff = { xzero, 0.0, 0.0 };
        a_ctr_grid = plt.addHexGridVisual (RD.hg, spatOff, zeromap, ctr_scaling);
        xzero += RD.hg->width();
    }

    if (plot_dr == true) {
        spatOff = { xzero, 0.0, 0.0 };
#if 0
        dr_grid = plt.addHexGridVisual (RD.hg, spatOff, zeromap, ctr_scaling);
#endif
        xzero += RD.hg->width();
    }

    // guidance expression
    if (plot_guide) {
        spatOff = { xzero, 0.0, 0.0 };
        float _m = 0.8; float _c = 0.0;
        const array<float, 4> guide_scaling = { 0.0f, 0.0f, _m, _c };
        // Plot gradients of the guidance effect g.
        for (unsigned int j = 0; j<RD.M; ++j) {
            plt.addHexGridVisual (RD.hg, spatOff, RD.rho[j], guide_scaling);
            spatOff[1] += RD.hg->depth();
        }
        xzero += RD.hg->width();
    }


    // Now plot fields and redraw display
    if (plot_guidegrad) {
        spatOff = { xzero, 0.0, 0.0 };
        for (unsigned int j = 0; j<RD.M; ++j) {

            // gradient of guidance expression
            vector<vector<FLT> > gx = separateVectorField (RD.g[j], 0);
            vector<vector<FLT> > gy = separateVectorField (RD.g[j], 1);
            FLT ming = 1e7;
            FLT maxg = -1e7;
            if (plot_guidegrad) {
                // Determine scale of gx and gy so that a common scale can be
                // applied to both gradient_x and gradient_y.
                for (unsigned int hi=0; hi<RD.nhex; ++hi) {
                    Hex* h = RD.hg->vhexen[hi];
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
            const array<float, 4> guidegrad_scaling = { 0.0f, 0.0f, gg_m, gg_c };

            // Create the grids
            plt.addHexGridVisual (RD.hg, spatOff, gx[j], guidegrad_scaling);
            spatOff[0] += RD.hg->width();
            plt.addHexGridVisual (RD.hg, spatOff, gy[j], guidegrad_scaling);
            spatOff[0] -= RD.hg->width();
            spatOff[1] += RD.hg->depth();
        }
        xzero += RD.hg->width() + 2.0f;
    }

    // Saving of t=0 images in log folder
    if ((RD.M > 0 && plot_guide) || plot_a) {
        savePngs (logpath, "sim", 0, plt);
    }

    // if using plotting, then set up the render clock
    steady_clock::time_point lastrender = steady_clock::now();

#endif // COMPILE_PLOTTING

    // Start the loop
    bool finished = false;
    while (finished == false) {
        // Step the model
        RD.step();

#ifdef COMPILE_PLOTTING

        if ((RD.stepCount % plotevery) == 0) {
            DBG2("Plot at step " << RD.stepCount);
            // Do a plot of the ctrs as found.
            vector<FLT> ctrmap = ShapeAnalysis<FLT>::get_contour_map (RD.hg, RD.c, RD.contour_threshold);

            if (plot_contours) {
                plt.updateHexGridVisual (c_ctr_grid, ctrmap, ctr_scaling);
            }

            if (plot_a_contours) {
                vector<FLT> actrmap = ShapeAnalysis<FLT>::get_contour_map (RD.hg, RD.a, RD.contour_threshold);
                plt.updateHexGridVisual (a_ctr_grid, actrmap, ctr_scaling);
            }

            if (plot_a) {
                for (unsigned int i = 0; i<RD.N; ++i) {
                    plt.updateHexGridVisual (agrids[i], RD.a[i], scaling);
                }
            }
            if (plot_c) {
                for (unsigned int i = 0; i<RD.N; ++i) {
                    plt.updateHexGridVisual (cgrids[i], RD.c[i], scaling);
                }
            }
            if (plot_n) {
                plt.updateHexGridVisual (ngrid, RD.n, scaling);
            }
#if 0
            if (plot_dr) {
                plt.updateHexGridVisual (dr_grid, RD.regions, ctr_scaling);
            }
#endif
            // Save to PNG
            if (vidframes) {
                savePngs (logpath, "sim", framecount, plt);
                ++framecount;
            } else {
                savePngs (logpath, "sim", RD.stepCount, plt);
            }
        }

        // rendering the graphics.
        steady_clock::duration sincerender = steady_clock::now() - lastrender;
        if (duration_cast<milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
            glfwPollEvents();
            plt.render();
            lastrender = steady_clock::now();
        }
#endif // COMPILE_PLOTTING

        // Save data every 'logevery' steps
        if (RD.stepCount == 1 || (RD.stepCount % logevery) == 0) {
            DBG ("Logging data at step " << RD.stepCount);
            RD.save();

            // Fixme. Save the hex contours in their own file. Each Hex has a save() method.
            vector<list<Hex> > sv_ctrs = ShapeAnalysis<FLT>::get_contours (RD.hg, RD.c, RD.contour_threshold);

            // If spatial analysis, then add line here to do it
        }

        if (RD.stepCount > steps) {
            finished = true;
        }
    }

    // Save out the sums.
    RD.savesums();

    // Before saving the json, we'll place any additional useful info
    // in there, such as the FLT. If float_width is 4, then
    // results were computed with single precision, if 8, then double
    // precision was used. Also save various parameters from the RD system.
    conf.set ("float_width", (unsigned int)sizeof(FLT));
    string tnow = morph::Tools::timeNow();
    conf.set ("sim_ran_at_time", tnow.substr(0,tnow.size()-1));
    conf.set ("hextohex_d", RD.hextohex_d);
    conf.set ("D", RD.get_D());
    conf.set ("k", RD.k);
    conf.set ("dt", RD.get_dt());
    // Call our function to place git information into root.
    //morph::Tools::insertGitInfo (conf.root, "sim/");
    conf.insertGitInfo ("sim/");
    // Store the binary name and command argument into root, too.
    if (argc > 0) { conf.set("argv0", argv[0]); }
    if (argc > 1) { conf.set("argv1", argv[1]); }

    // We'll save a copy of the parameters for the simulation in the log directory as params.json
    const string paramsCopy = logpath + "/params.json";
    conf.write (paramsCopy);
    if (conf.ready == false) {
        cerr << "Warning: Something went wrong writing a copy of the params.json: " << conf.emsg << endl;
    }

    // Extract contours
    vector<list<Hex> > ctrs = ShapeAnalysis<FLT>::get_contours (RD.hg, RD.c, RD.contour_threshold);
    {
        // Write each contour to a contours.h5 file
        stringstream ctrname;
        ctrname << logpath << "/contours.h5";
        HdfData ctrdata(ctrname.str());
        unsigned int nctrs = ctrs.size();
        ctrdata.add_val ("/num_contours", nctrs);
        for (unsigned int ci = 0; ci < nctrs; ++ci) {
            vector<FLT> vx, vy;
            auto hi = ctrs[ci].begin();
            while (hi != ctrs[ci].end()) {
                vx.push_back (hi->x);
                vy.push_back (hi->y);
                ++hi;
            }
            stringstream ciss;
            ciss << ci;
            string pth = "/x" + ciss.str();
            ctrdata.add_contained_vals (pth.c_str(), vx);
            pth[1] = 'y';
            ctrdata.add_contained_vals (pth.c_str(), vy);

            // Generate hex grids from contours to obtain the size of the region enclosed by the contour
            HexGrid* hg1 = new HexGrid (RD.hextohex_d, RD.hexspan, 0, morph::HexDomainShape::Boundary);
            hg1->setBoundary (ctrs[ci]);
            pth[1] = 'n';
            ctrdata.add_val(pth.c_str(), hg1->num());
            delete hg1;
        }

        // Also extract the boundary of the main, enclosing hexgrid and write that.
        list<Hex> outerBoundary = RD.hg->getBoundary();
        vector<FLT> vx, vy;
        auto bi = outerBoundary.begin();
        while (bi != outerBoundary.end()) {
            vx.push_back (bi->x);
            vy.push_back (bi->y);
            ++bi;
        }
        ctrdata.add_contained_vals ("/xb", vx);
        ctrdata.add_contained_vals ("/yb", vy);
    }

#ifdef COMPILE_PLOTTING
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    plt.keepOpen();
#endif

    return 0;
};
