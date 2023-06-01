# Optimisation topologique d’un rotor de MRV pour maximiser le couple moyen : un problème mal posé
Codes associés à l'article "Optimisation topologique d’un rotor de MRV pour maximiser le couple moyen : un problème mal posé", SGE 2023, Cherriere T., Gauthey T., Hlioui S., Hage-Hassan M.

![GitHub license](https://img.shields.io/github/license/tcherrie/SGE23_Cherriere_Gauthey_Hlioui_HageHassan)
[![DOI](https://zenodo.org/badge/642540718.svg)](https://zenodo.org/badge/latestdoi/642540718)


## 1) Contenu : 

- `README.RD` : ce fichier
- `module_optim_topo.py` : module python contenant dertaines macros
- `Lamine lineaire.ipynb` : notebook sur les propriété d'un laminé d'ordre 1 homogénéisé linéaire
- `Lamine non-lineaire.ipynb` : notebook sur les propriété d'un laminé d'ordre 1 homogénéisé non-linéaire
- `Suite minimisante lineaire.ipynb`: notebook qui trace l'évolution de la compliance en fonction de la finesse des couches de laminé (linéaire)
- `Suite minimisante non-lineaire.ipynb`: notebook qui trace l'évolution de la compliance en fonction de la finesse des couches de laminé (non-linéaire)
- `Topopt_lineaire.ipynb` : notebook qui propose une optimisation topologique d'un rotor de MRV (fer considéré linéaire)
- `Topopt_non_lineaire.ipynb` : notebook qui propose une optimisation topologique d'un rotor de MRV (fer considéré non-linéaire)

## 2) Instructions d'installation locale

1. Télécharger et décompresser le contenu de ce [répertoire github](https://github.com/tcherrie/SGE23_Cherriere_Gauthey_Hlioui_HageHassan) **sur votre bureau**. Pour cela cliquer sur "Code", puis "Download ZIP".
2. Télécharger et installer [Miniconda](https://docs.conda.io/en/latest/miniconda.html) avec les options par défaut correspondant à votre système
3. Ouvrir une console Anaconda Prompt (sur Windows, taper dans la barre de recherche en bas à gauche "Anaconda prompt" ; sur Linux et Mac on utilisera directement le terminal)

Désormais, on écrira la suite d'instructions suivante sur la console (vous pouvez copier les instructions et les coller dans la console par un clic droit ; appuyer sur Entrée après avoir écrit l'instruction pour l'exécuter ; appuyer sur "y" pour confirmer lorsque cela vous sera demandé) : 

4. Créer un nouvel environnement "topopt_SGE23" : `conda create -n topopt_SGE23 python=3`
5. Activer l'environnement : `conda activate topopt_SGE23`
6. Installer des packages nécessaires (500Mb) : `conda install jupyter numpy matplotlib scipy` (appuyer sur y + entrée pour confirmer l'installation). Si un blocage apparait, fermer la console et recommencer les étapes 3, 5 et 6. Les packages déjà téléchargés ne seront pas retéléchargés, ce qui fluidifiera le processus.
7. Installer NGSolve (300Mb) : `pip install ngsolve`
8. Installer des extensions de visualisation sur les notebooks :
 - `pip install --upgrade webgui_jupyter_widgets`
 - `jupyter nbextension install --user --py webgui_jupyter_widgets`
 - `jupyter nbextension enable --user --py webgui_jupyter_widgets`

**Si une erreur apparaît à l'étape 8**, essayer `jupyter nbextension install webgui_jupyter_widgets` ; sinon passer directement à l'étape 9.

9. Lancer Jupyter : `jupyter notebook`
10. Dans l'explorateur de fichiers qui s'ouvre, sélectionner "Desktop", puis "SGE23_Cherriere_Gauthey_Hlioui_HageHassan"

Vous pouvez désormais exécuter et interagir avec les notebooks.

## 3) License

Copyright (C) 2023 Théodore CHERRIERE (theodore.cherriere@ens-paris-saclay.fr), Thomas GAUTHEY (thomas.gauthey@centralesupelec.fr), Sami HLIOUI (sami.hlioui@cyu.fr), Maya HAGE-HASSAN (maya.hage-hassan@centralesupelec.fr)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
