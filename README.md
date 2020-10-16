# RetinoTectal
A model of the self-organization of the retino-tectal projection

To build:

```bash
# If you didn't already, then clone RetinoTectal:
git clone https://github.com/ABRG-Models/RetinoTectal.git
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
