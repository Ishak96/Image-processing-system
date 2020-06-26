# Système complet de suivi d’objets

Réalisation d'un système complet de suivi d'objet dans une séquence d'images couleur, en utilisant la méthode de suivi par sélection de vignettes autour des points d'intérêt ainsi que la méthode du suivi Milti-critères.

## Pour commencer

Récupérer le répértoire avec :

```bash
$ git clone https://github.com/Ishak96/image-processing-system.git
```

Ou bien directement du site.

### Dans le dossier videos

Vous allez trouver des séquences vidéo enregistrées illustrant les différents résultats obtenus sur les deux séquences d'images, ainsi qu'un exemple complet d'utilisation de l'application.

Au niveau de l'architecture du code, j'ai découpé en plusieurs dossiers.

- include/ comprenant les entêtes
- src/ pour les fichiers .c
- bin/ pour les exécutables
- test/ comprenant les tests des différentes parties

### Pré-requis

- La bibliothèque NRC
- GTK+ 3

### Installation

Pour que l'application fonctionne sans problèmes de compilation il faut que vous installiez la bibliothèque GTK+ 3.

Nous allons commencer par mettre vos paquets à jour, ouvrez donc la console et tapez :
```bash
$ sudo apt-get update
```

Maintenant, installez le paquet de développement (qui contient les bibliothèques).
C'est toujours aussi simple, il suffit de faire :
```bash
$ sudo apt-get install libgtk2.0-dev
```

Vous devriez normalement avoir tout ce qu'il faut pour développer avec GTK, mais il faut s'assurer que vous disposez bien du Runtime :
```bash
$ sudo apt-get install libgtk2.0-0
```

## Démarrage

Déplacez vous vers le dossier NRC :
```bash
$ cd NRC/
```

Lancez le makefile pour qu'il génère les exécutable des tests ainsi que le main :
```bash
$ make && make test
```

Pour lancer un test :
```bash
$ bin/nom_test
```

Pour lancer l'application :
```bash
$ bin/object_tracking "LABEL"
```

Le résultat sera stocké dans le dossier /Images/results/object_tracking/

Précisez Le "LABEL" -> Lbox ou Fomd

Pour lancer que le suivi par sélection de vignette autour des points d'intérêt :

Pour lancer l'application :
```bash
$ bin/vignette_suivi_objet
```

Le résultat sera stocké dans le dossier /Images/results/morphoMath/extraxtion/movementDetectiontRefImage/Lbox/

## Fabriqué avec

_exemples :_
* [GTK+ 3](https://developer.gnome.org/gtk3/stable/) - The GIMP Toolkit, anciennement GTK+
* [NRC] - bibliothèque pour le traitement d'image

## Auteur
* **Ishak Ayad** [@Ishak96](https://github.com/Ishak96)
