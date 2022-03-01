/*
 * Script-like program to read data from annealling runs and visualise the results
 */

#include <sstream>
#include <string>
#include <stdexcept>
#include <map>
#include <morph/vVector.h>
#include <morph/Config.h>
#include <morph/HdfData.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/ScatterVisual.h>
#include <morph/TriaxesVisual.h>

// Dimensionality.
#ifndef D
# define D 3
#endif

int main (int argc, char** argv)
{
    // Load data, based on command line arg
    if (argc < 2) { throw std::runtime_error ("Provide path to data file"); }

    // Filename something like: path/anneal1_ee_GJ_20211015_145811.h5
    std::string fname(argv[1]);
    morph::Tools::stripUnixPath (fname);
    std::string::size_type strt = fname.find ("anneal1_") + 8;
    std::string::size_type nd = fname.find ("_20");
    std::string modelid = fname.substr (strt, nd-strt);
    std::string model = std::string("configs/a1/m_") + modelid + std::string(".json");;
    std::cout << "model is " << model << std::endl;

    // Open JSON file (for later saving)
    morph::Config mdl(model);

    // Read the HDF5 into data structures
    morph::HdfData data(argv[1], morph::FileAccess::ReadOnly);
    morph::vVector<morph::Vector<float, D>> param_hist_accepted;
    // To hold param_hist_accepted parameters in VisualModel coordinates for plotting
    morph::vVector<morph::Vector<float, D>> param_hist_accepted_vm_coords;
    morph::vVector<float> f_param_hist_accepted;
    morph::vVector<morph::Vector<float, D>> param_hist_rejected;
    morph::vVector<morph::Vector<float, D>> param_hist_rejected_vm_coords;
    morph::vVector<float> f_param_hist_rejected;
    morph::vVector<float> T_k_hist;
    morph::vVector<float> T_cost_hist;
    morph::vVector<float> f_x_hist;
    morph::vVector<float> f_x_best_hist;
    morph::vVector<float> x_best;
    data.read_contained_vals ("/param_hist_accepted", param_hist_accepted);
    data.read_contained_vals ("/param_hist_accepted", param_hist_accepted_vm_coords);
    data.read_contained_vals ("/f_param_hist_accepted", f_param_hist_accepted);
    data.read_contained_vals ("/param_hist_rejected", param_hist_rejected);
    data.read_contained_vals ("/param_hist_rejected", param_hist_rejected_vm_coords);
    data.read_contained_vals ("/f_param_hist_rejected", f_param_hist_rejected);
    data.read_contained_vals ("/T_k_hist", T_k_hist);
    data.read_contained_vals ("/T_cost_hist", T_cost_hist);
    data.read_contained_vals ("/f_x_hist", f_x_hist);
    data.read_contained_vals ("/f_x_best_hist", f_x_best_hist);
    data.read_contained_vals ("/x_best", x_best);

    std::cout << "Total parameter sets examined: "
              << (param_hist_accepted.size() + param_hist_rejected.size()) << std::endl;

    morph::Vector<float, D> range_min;
    morph::Vector<float, D> range_max;
    std::vector<std::string> pnames;
    data.read_contained_vals ("/range_max", range_max);
    data.read_contained_vals ("/range_min", range_min);
    for (size_t i = 1; i <= D; ++i) {
        std::string pn = std::string("/param_name_") + std::to_string(i);
        std::string pname("");
        data.read_string (pn.c_str(), pname);
        std::cout << "parameter name: " << pname << std::endl;
        pnames.push_back (pname);
    }

    morph::Vector<float, D> range_diff = range_max - range_min;
    morph::Vector<float, D> range_offs = range_min / range_diff;
    // This converts param_hist_rejected to model coords
    for (auto& phr : param_hist_rejected_vm_coords) {
        phr /= range_diff;
        phr -= range_offs;
    }
    for (auto& pha : param_hist_accepted_vm_coords) {
        //pha += range_min;
        pha /= range_diff;
        pha -= range_offs;
    }

    // Use a map to order all rejected and accepted parameter sets
    std::map<float, morph::Vector<float, D>> mapped_params; // Parameter space
    std::map<float, morph::Vector<float, D>> mapped_params_vm_coords; // Scaled parameter space
    for (size_t i = 0; i < param_hist_accepted.size(); ++i) {
        mapped_params[f_param_hist_accepted[i]] = param_hist_accepted[i];
        mapped_params_vm_coords[f_param_hist_accepted[i]] = param_hist_accepted_vm_coords[i];
    }
    for (size_t i = 0; i < param_hist_rejected.size(); ++i) {
        mapped_params[f_param_hist_rejected[i]] = param_hist_rejected[i];
        mapped_params_vm_coords[f_param_hist_rejected[i]] = param_hist_rejected_vm_coords[i];
    }

    // Now, for each of the n best param sets, find the worse parameter set within a
    // certain radius *in scaled space*. Record the
    float r = 0.1f;
    std::map<float, morph::Vector<float, D>>::const_iterator pi = mapped_params_vm_coords.begin();
    size_t j = 0;
    // container for the 'biggest within range', keyed by the objective of the *good* parameters
    std::map<float, std::pair<float, morph::Vector<float, D>>> biggest_nearby;

    while (pi != mapped_params_vm_coords.end() && j++ < 5) {
        morph::Vector<float, D> smallparam = pi->second;
        float smallobj = pi->first;
        // Run though the params finding the largest one within radius
        float maxg = std::numeric_limits<float>::lowest();
        morph::Vector<float, D> maxgparams;
        float maxgobj = 0.0f;
        maxgparams.zero();
        for (size_t i = 0; i < param_hist_accepted_vm_coords.size(); ++i) {
            morph::Vector<float, D> testparam = param_hist_accepted_vm_coords[i];
            float testobj = f_param_hist_accepted[i];
            float d = (testparam-smallparam).length();
            if (d < r) {
                // Close enough to consider. What's the gradient?
                float g = std::abs(testobj - smallobj) / d;
                maxgparams = g > maxg ? testparam : maxgparams;
                maxgobj = g > maxg ? testobj : maxgobj;
                maxg = g > maxg ? g : maxg;
            }
        }
        for (size_t i = 0; i < param_hist_rejected_vm_coords.size(); ++i) {
            morph::Vector<float, D> testparam = param_hist_rejected_vm_coords[i];
            float testobj = f_param_hist_rejected[i];
            float d = (testparam-smallparam).length();
            if (d < r) {
                // Close enough to consider. What's the gradient?
                float g = std::abs(testobj - smallobj) / d;
                maxgparams = g > maxg ? testparam : maxgparams;
                maxgobj = g > maxg ? testobj : maxgobj;
                maxg = g > maxg ? g : maxg;
            }
        }
        // Write maxgparams into container using smallobj as key
        biggest_nearby[smallobj] = std::make_pair(maxgobj, maxgparams);
        ++pi;
    }

    // How about I write out JSON config files to run the agent for the best parameter set(s)?
    std::cout << "Best params: " << x_best << std::endl;

    auto mp = mapped_params.begin(); // Should be the best
    mp++; // Next best
    mdl.set (pnames[0], mp->second[0]);
    mdl.set (pnames[1], mp->second[1]);
    mdl.set (pnames[2], mp->second[2]);
    mdl.write (modelid+std::string("_x_2nd.json"));
    mp++; // Next best
    mdl.set (pnames[0], mp->second[0]);
    mdl.set (pnames[1], mp->second[1]);
    mdl.set (pnames[2], mp->second[2]);
    mdl.write (modelid+std::string("_x_3rd.json"));

    morph::vVector<float> simtime (T_k_hist.size());
    simtime.linspace (0.0, (float)T_k_hist.size());

    // Set up the visualisation
    morph::Visual v (1920, 1080, "Optimization progress");
    v.zNear = 0.001;
    v.setSceneTransZ (-3.0f);
    v.lightingEffects (true);
    morph::Vector<float, 3> offset = { -0.7, 0.0, 0.0 };

    // First a scatter plot that can be updated. Just using a ScatterVisual for this.
    morph::ScatterVisual<float>* sv = new morph::ScatterVisual<float> (v.shaderprog, offset);
    sv->radiusFixed = 0.002f;
    sv->colourScale.compute_autoscale (0, 30);
    sv->cm.setType (morph::ColourMapType::Jet);
    sv->setDataCoords ((std::vector<morph::Vector<float, 3>>*)&param_hist_accepted_vm_coords);
    sv->setScalarData ((std::vector<float>*)&f_param_hist_accepted);
    sv->sizeFactor = 0.1f; //0.05f is 1.0f/300.0f;
    sv->finalize();
    v.addVisualModel (sv);

    morph::ScatterVisual<float>* sv2 = new morph::ScatterVisual<float> (v.shaderprog, offset);
    sv2->radiusFixed = 0.002f;
    sv2->colourScale.compute_autoscale (0, 30);
    sv2->cm.setType (morph::ColourMapType::Plasma);
    sv2->setDataCoords ((std::vector<morph::Vector<float, 3>>*)&param_hist_rejected_vm_coords);

    static constexpr float scatter_max_sz = 1.0f; // How big should the bad blobs be allowed to get?
    morph::vVector<float> log_fphr = f_param_hist_rejected.log();
    sv2->setScalarData ((std::vector<float>*)&log_fphr);
    sv2->sizeFactor = sv->sizeFactor;
    sv2->finalize();
    v.addVisualModel (sv2);

    // Show the best one or best few:
    morph::ScatterVisual<float>* sv3 = new morph::ScatterVisual<float> (v.shaderprog, offset);
    sv3->radiusFixed = 0.01f;
    sv3->colourScale.compute_autoscale (0, 1);
    sv3->cm.setType (morph::ColourMapType::Jet);
    sv3->finalize();
    // Show the 5 best and their neighbouring biggest gradient partners
    j = 1;
    for (auto bn : biggest_nearby) {
        // bn.first = obj/key, bn.second.second is partner bn.second.first is partner objective
        morph::Vector<float, D> coord = mapped_params_vm_coords[bn.first];
        std::cout << "Best objective: " << bn.first << ", neighbour obj: " << bn.second.first << std::endl;
        sv3->add (coord, 1.0f, bn.first/100.0f); // one of the best
        sv3->add (bn.second.second, 0.0f, bn.second.first/1000.0f); // the most changed params nearby

        // Also write out jsons. One for the 'best' value
        coord = (coord + range_offs) * range_diff;
        mdl.set (pnames[0], coord[0]);
        mdl.set (pnames[1], coord[1]);
        mdl.set (pnames[2], coord[2]);
        mdl.set ("fobj", bn.first);
        mdl.write (std::string("./log/anneal1/m_")+modelid+("_best_")+std::to_string(j)+(".json"));
#if 0
        // And one for the params that are near to the 'best' value and have the greatest change/distance from it
        morph::Vector<float, D> coord2 = (bn.second.second + range_diff) * range_max;
        mdl.set (pnames[0], coord2[0]);
        mdl.set (pnames[1], coord2[1]);
        mdl.set (pnames[2], coord2[2]);
        mdl.set ("fobj", bn.second.first);
        mdl.write (std::string("./log/anneal1/m_")+modelid+("_nearby_")+std::to_string(j)+(".json"));
#endif
        ++j;
    }

    v.addVisualModel (sv3);

    morph::TriaxesVisual<float>* tav = new morph::TriaxesVisual<float> (v.shaderprog, v.tshaderprog, offset);
    tav->axisstyle = morph::axisstyle::L;
    tav->input_min = range_min;
    tav->input_max = range_max;
    tav->xlabel = pnames[0];
    tav->ylabel = pnames[1];
    tav->zlabel = pnames[2];
    tav->finalize();
    v.addVisualModel (tav);

    offset[0] += 2.0f;
    // Add a graph to track T_i and T_cost
    morph::GraphVisual<float>* graph1 = new morph::GraphVisual<float> (v.shaderprog, v.tshaderprog, offset);
    graph1->twodimensional = false;
    graph1->policy = morph::stylepolicy::lines;
    graph1->ylabel = "log(T)";
    graph1->xlabel = "Anneal time";
    graph1->setdata (simtime, T_k_hist, "Tparam");
    graph1->setdata (simtime, T_cost_hist, "Tcost");
    graph1->finalize();
    v.addVisualModel (graph1);

    offset[0] += 1.4f;
    morph::GraphVisual<float>* graph2 = new morph::GraphVisual<float> (v.shaderprog, v.tshaderprog, offset);
    graph2->twodimensional = false;
    graph2->policy = morph::stylepolicy::lines;
    graph2->ylabel = "obj value";
    graph2->xlabel = "Anneal time";
    graph2->setdata (simtime, f_x_hist, "f_x");
    graph2->setdata (simtime, f_x_best_hist, "f_x_best");
    graph2->finalize();
    v.addVisualModel (graph2);
#if 0
    offset[0] += 1.4f;
    morph::GraphVisual<double>* graph3 = new morph::GraphVisual<double> (v.shaderprog, v.tshaderprog, offset);
    graph3->twodimensional = false;
    graph3->setlimits (0, 1000, -1.0f, 100.0f);
    graph3->policy = morph::stylepolicy::lines;
    graph3->ylabel = "obj value";
    graph3->xlabel = "Anneal time";
    graph3->prepdata ("f_x_cand");
    graph3->finalize();
    v.addVisualModel (graph3);

    // Text labels to show additional information that might update
    morph::Vector<float> lpos = {-0.08f, 0.03f, 0.0f};
    morph::VisualTextModel* fps_tm;
    v.addLabel ("Unset", lpos, fps_tm);

    // Fixed text labels
    std::stringstream ss;
    lpos[1] -= 0.02f;
    ss << "reanneal_after_steps = " << optimiser->reanneal_after_steps;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "temperature_ratio_scale = " << optimiser->temperature_ratio_scale;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "temperature_anneal_scale = " << optimiser->temperature_anneal_scale;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "cost_parameter_scale_ratio = " << optimiser->cost_parameter_scale_ratio;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "acc_gen_reanneal_ratio = " << optimiser->acc_gen_reanneal_ratio;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "objective_repeat_precision = " << optimiser->objective_repeat_precision;
    v.addLabel (ss.str(), lpos);

    lpos[1] -= 0.02f; ss.str("");
    ss << "f_x_best_repeat_max = " << optimiser->f_x_best_repeat_max;
    v.addLabel (ss.str(), lpos);
#endif

    v.keepOpen();

    return 0;
}
