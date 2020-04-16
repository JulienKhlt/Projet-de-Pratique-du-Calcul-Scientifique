# Projet de Pratique du Calcul Scientifique

Modélisation d'une épidémie sous Julia

Actuellement il y a deux fichiers :

- SIR_models.jl
- Diffusion_model.jl

# Table of Contents

- [Projet de Pratique du Calcul Scientifique](#projet-de-pratique-du-calcul-scientifique)
- [Table of Contents](#table-of-contents)
  - [Modèle SIR](#mod%c3%a8le-sir)
    - [TODO](#todo)
  - [Modèle diffusion](#mod%c3%a8le-diffusion)
    - [TODO](#todo-1)
  - [Développement](#d%c3%a9veloppement)
    - [TODO](#todo-2)
  - [Requirement](#requirement)


## Modèle SIR

*Premier modèle simple qui se base sur un système d'équation relativement basique.* <br>
On note **S** la population susceptible d'avoir le virus, <br>
**I** la population infectée par le virus, <br>
et **R** la population guérit ou morte du virus, <br> 
**On a ainsi les équations suivantes:**
$$
    \frac{dS}{dt} = - \frac{\beta SI}{N}
$$
$$
    \frac{dI}{dt} = \frac{\beta SI}{N} - \gamma I
$$
$$
    \frac{dR}{dt} = \gamma I
$$

avec $N = S + I + R$ <br>
**Ceci correspond au modèle SIR.** <br>

*Lorsqu'executé le document créera des courbes grace à PyPlot qu'il faut avoir installé.*

### TODO

- [ ] Existence et Unicité de la solution
- [ ] Point fixes
- [ ] Stabilité
- [ ] Invariants
- [ ] Comportement global.

## Modèle diffusion

On reprend les mêmes équations que précédement mais en y ajoutant de la diffusion.
*Lorsqu'executé le document créera des courbes grace à PyPlot qu'il faut avoir installé.* <br>
**Attention il faut que ![\frac{Ddt}{dx^2} < \frac{1}{2}](https://render.githubusercontent.com/render/math?math=%5Cfrac%7BDdt%7D%7Bdx%5E2%7D%20%3C%20%5Cfrac%7B1%7D%7B2%7D)**

### TODO

- [ ] Conditions limites
- [ ] Dérivation de cette ́equation
- [ ] Rôle de D
- [ ] Points fixes
- [ ] Invariants
- [ ] Comportement global des solutions
- [ ] Comment valider votre résultat ?
- [ ] Faire une  ́etude de stabilité de la méthode d’Euler explicite sur $\dot{x}=−\lambda x,λ >0$
- [ ] Puis sur ̇$\dot{x}=−Ax$ *où A est une matrice symétrique définie positive.* 
- [ ] En déduire une relation entre $∆t$ et $∆x$ pour le cas de l’équation de la chaleur. Que se passe-t-il si cette condition n’est pas vérifiée ? 
- [ ] A quelle erreur  vous  attendez-vous  en  fonction  de  $∆x$ et  $∆t$? Vérifier numériquement. 
- [ ] Implémentez la méthode d’Euler implicite pour cette  ́equation, et résolvez le système obtenir par la méthode de Newton.  
- [ ] Peut-on garantir la convergence de la méthode de Newton?  
- [ ] Quelle est la complexité de cet algorithme par rapport au nombre de points $N$?  Utiliser une méthode d’inversion itérative. 
- [ ] Quel est le coût par itération en fonction de $N$?  Quel est le coût total ?

## Développement

### TODO
- [ ] Maths: Etude d’existence et d’unicité de solutions (recherche bibliographique),analyse de stabilité d’EDP...
- [ ] Méthodes d’intégration: essayer d’autres méthodes (ordre supérieur, méthodes avec contrôle d’erreur...), utiliser des bibliothèques existantes...
- [ ] Résolution d’équations:  utiliser une méthode de Krylov, utiliser des bibliothèques existantes...
- [ ] Modélisation: 
   - Changer de classe de modèle (modèle stochastique, à agents...),
   - Complexifier celui-ci (SEIR), 
   - Utiliser un modèle pour illustrer ou répondre à une question pratique (comment confiner une population,
   - Quelle stratégie pour minimiser le nombre de morts, impact ́economique...)
- [ ] Données:   calibrer  le  modèle,  tenter  de  prédire  l’évolution, à partir  de quand dans une épidémie peut-on raisonnablement extrapoler une courbe...
- [ ] Informatique: optimiser le code, le paralléliser, faire une interface graphique..


## Requirement

- PyPlot
