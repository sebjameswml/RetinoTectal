# RetinoTectal
A model of the self-organization of the retino-tectal projection

To build, first install morphologica dependences, as described in morphologica's [README.install.linux.md](https://github.com/ABRG-Models/morphologica/blob/master/README.install.linux.md) or [README.install.mac.md](https://github.com/ABRG-Models/morphologica/blob/master/README.install.mac.md), the follow this recipe:

```bash
# If you didn't already, then clone RetinoTectal:
git clone https://github.com/sebjameswml/RetinoTectal.git
cd RetinoTectal
# Clone the header-only morphologica in-tree:
git clone https://github.com/ABRG-Models/morphologica.git
mkdir build
pushd build
cmake ..
make -j4 # or however many cores you have
# (no need to install, you'll run the simulations in place)
popd
```

There are actually two simulation codes in this repository. The
original work was a reaction-diffusion style model found in
sim/rd. The final model was an agent-based model in sim/agent/. The
agent-based model was the one used as the basis for the paper.

To run the simulations described in the paper, you can run:

./scripts/paper_figures.sh
