#include <sstream>
#include <string>
#include <stdexcept>
#include <morph/vVector.h>
#include <morph/HdfData.h>
#include <morph/Visual.h>
#include <morph/GraphVisual.h>
#include <morph/ScatterVisual.h>
#include <morph/TriaxesVisual.h>

int main (int argc, char** argv)
{
    // Load data, based on command line arg
    if (argc < 2) { throw std::runtime_error ("Provide path to data file"); }

    morph::HdfData data(argv[1], morph::FileAccess::ReadOnly);
    morph::vVector<morph::Vector<float, 3>> param_hist_accepted;
    morph::vVector<float> f_param_hist_accepted;
    morph::vVector<morph::Vector<float, 3>> param_hist_rejected;
    morph::vVector<float> f_param_hist_rejected;
    morph::vVector<float> T_k_hist;
    morph::vVector<float> T_cost_hist;
    morph::vVector<float> f_x_hist;
    morph::vVector<float> f_x_best_hist;
    data.read_contained_vals ("/param_hist_accepted", param_hist_accepted);
    data.read_contained_vals ("/f_param_hist_accepted", f_param_hist_accepted);
    data.read_contained_vals ("/param_hist_rejected", param_hist_rejected);
    data.read_contained_vals ("/f_param_hist_rejected", f_param_hist_rejected);
    data.read_contained_vals ("/T_k_hist", T_k_hist);
    data.read_contained_vals ("/T_cost_hist", T_cost_hist);
    data.read_contained_vals ("/f_x_hist", f_x_hist);
    data.read_contained_vals ("/f_x_best_hist", f_x_best_hist);

    // This code section just hacked in until I regenerate data with range_max/range_min in the code
    morph::vVector<morph::Vector<float,2>> param_ranges;
    morph::Vector<float, 3> range_min;
    morph::Vector<float, 3> range_max;
    std::vector<std::string> pnames;
    for (size_t i = 1; i <= 3; ++i) {
        std::string pn = std::string("/param_name_") + std::to_string(i);
        std::string pname("");
        data.read_string (pn.c_str(), pname);
        std::cout << "parameter name: " << pname << std::endl;
        pnames.push_back (pname);
        if (pname == "r_c") {
            param_ranges.push_back ({0.001f, 0.5f});
            range_min[i-1] =  (0.001f);
            range_max[i-1] =  (0.5f);
        } else if (pname == "r_i") {
            param_ranges.push_back ({0.001f, 0.5f});
            range_min[i-1] =  (0.001f);
            range_max[i-1] =  (0.5f);
        } else if (pname == "r_j") {
            param_ranges.push_back ({0.001f, 0.5f});
            range_min[i-1] =  (0.001f);
            range_max[i-1] =  (0.5f);
        } else if (pname == "s") {
            param_ranges.push_back ({0.01f, 0.99f});
            range_min[i-1] =  (0.01f);
            range_max[i-1] =  (0.99f);
        } else if (pname == "m_g") {
            param_ranges.push_back ({0.0001f, 0.01f});
            range_min[i-1] =  (0.0001f);
            range_max[i-1] =  (0.01f);
        } else if (pname == "m_c") {
            param_ranges.push_back ({0.01f, 0.1f});
            range_min[i-1] =  (0.01f);
            range_max[i-1] =  (0.1f);
        } else if (pname == "m_i") {
            param_ranges.push_back ({0.01f, 0.8f});
            range_min[i-1] =  (0.01f);
            range_max[i-1] =  (0.8f);
        } else if (pname == "m_j") {
            param_ranges.push_back ({0.00001f, 0.001f});
            range_min[i-1] =  (0.00001f);
            range_max[i-1] =  (0.001f);
        }
    }

    for (auto& phr : param_hist_rejected) { phr /= range_max; }
    for (auto& pha : param_hist_accepted) { pha /= range_max; }
    //param_hist_rejected /= range_max;

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
    sv->setDataCoords ((std::vector<morph::Vector<float, 3>>*)&param_hist_accepted);
    sv->setScalarData ((std::vector<float>*)&f_param_hist_accepted);
    sv->sizeFactor = 0.05f;//1.0f/300.0f;
    sv->finalize();
    v.addVisualModel (sv);

    morph::ScatterVisual<float>* sv2 = new morph::ScatterVisual<float> (v.shaderprog, offset);
    sv2->radiusFixed = 0.002f;
    sv2->colourScale.compute_autoscale (0, 30);
    sv2->cm.setType (morph::ColourMapType::Plasma);
    sv2->setDataCoords ((std::vector<morph::Vector<float, 3>>*)&param_hist_rejected);
    sv2->setScalarData ((std::vector<float>*)&f_param_hist_rejected);
    sv2->sizeFactor = sv->sizeFactor;
    sv2->finalize();
    v.addVisualModel (sv2);

    morph::TriaxesVisual<float>* tav = new morph::TriaxesVisual<float> (v.shaderprog, v.tshaderprog, offset);
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
