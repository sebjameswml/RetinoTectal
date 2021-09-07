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

whereas a model with chemo, and interaction woudl be

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
