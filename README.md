# BE Diffusion advection 2D

[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)

Programme Fortran de résolution numérique d’une équation aux dérivées partielles.
d’advection-diffusion par une méthode de type Volumes Finis.

Créé le 14/05/21 par Gabriel Darremont et Sam Alonso-Virissel.

## Pré-requis

- Système linux (Pour Windows, effacer les lignes 41 à 44 de _VTSWriter.f90_)
- GNU Fortran 7.5.0 ou supérieur (compilation)
- Paraview (visualisation)

## Installation

1. Extraire les sources dans un dossier et se placer avec un terminal dedans.
2. Exécuter ``make`` pour compiler les sources.

## Exécution

1. Ajuster les paramètres d'entrée dans _donnee.dat_.
2. Lancer le calcul avec la commande ``./prog.exe``.
3. Visualiser le résultat dans Paraview comme suit :
   1. ``File > Load State`` et ouvrir _sol.pvsm_.
   2. Visualiser l'évolution en fonction du temps avec ``Play``.

Le programme génère un fichier _sol.pvd_ à la racine ainsi qu'un maximum de 100 fichiers sol******.vts dans le dossier vts_file.
