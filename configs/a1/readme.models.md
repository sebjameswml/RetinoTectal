# Models

Model configs are prefixed 'm_'

'm_ee' models are 'exponential-exponential' models, referring to the
idea that both retinal and tectal expression patterns are of
exponential (and not linear or logarithmic form).

'm_el' models are 'exonential-logarithmic'. Here, the tectum has a
logarithm form in one or more ligand types. This is necessary where
the receptor-ligand interaction is attractive, rather than repulsive.

A model with 'axonaxon' in its name is one in which the axon-axon
receptor-ligand interaction provides the mechanism of competition for
space, rather than a distance based competition of unknown mechanism.

A model with 'chemoonly' is one in which competitions and axon-axon
interactions are switched off, leaving only the chemoaffinity
component.

A better scheme would be to mirror the letters used in the paper. So a model with chemoaffinity only might be

m_ee_G.json

whereas a model with chemo, and interaction would be

m_ee_GI.json

And a model with chemo and two types of axon-axon interaction:

m_ee_GIJ.json

where I've just decided that the second kind of interaction should be
called J, but that J should be the receptor-ligand interaction and I
should remain the receptor-receptor interaction that seems to work to
reproduce Reber/Brown results.

If I'm using the 'alternative' exponential curves, then the models could be:

m_eE_GIJ.json

or

m_eE_GCIJ.json


m_eel*.json files are exponential-exponential with one logarithmic
ligand on the tectum (modelling EphB).

m_lin* models are linear models - receptors and ligands vary linearly

m_linattr* models are linear models with attractive interactions.

If the m_XX tag ends with '2' then it;s a model with only 2 receptor-ligands per sheet instead of 4.

m_ql models are 'quadratic-linear'

m_random.json is a special model to randomly locate branches for performing statistical analysis

m_stoc.json is a special stochastic model, yet to be fully developed

m_geb.json sets up the Gebhardt et al model

## Autogeneration of models

Because there are so many permutations, I decided to write
make_models.py to generate all the json configs. That means if I
change a common parameter (like branches per axon, bpa) then I can
easily re-generate ALL of these files.
