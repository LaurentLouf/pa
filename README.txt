Configuration et execution testé exclusivement sous linux :

installer CMake 2.8+ au besoin
installer g++ au besoin
installer libgomp1 si cmake si nécessaire (librairie OpenMP)

dans le dossier ou se trouve ce README, tapez
cmake .
make
./ballDetection

Si un problème survient à propos de la librairie ocl et produit une erreur de segmentation, c'est que votre pilote de carte graphique n'est pas supporté, désactivez le define OCL dans imageProcessing.h

pour activer ou desactiver l'affichage faites de même avec le define DISPLAY

Si des problèmes apparaissent avec OpenMP vous pouvez le désactiver en commentant le define OMP dans lib.h

Il est possible de changer le fichier de capture, pour celà changer tout simplement le nom du fichier dans ballDetection.cpp à la ligne 93, si c'est une image activez le define PICTURE
