# Fonction ACP permet d'appliquer la méthode de l'analyse en composantes principales normée
#et non normée (indication par scale qui prend True par défaut)
#Elle prend en argument un tableau de données contenant variables quantitatives et qualitatives 
#puis garde seulement les quantitatives pour les calculs 
#Pour un tableau de données avec des valeurs manquantes, il y a un arret du programme 
#La fonction renvoi les objets suivants: 
#                     - Tableau des valeurs centrées (réduites)
#                     - matrice de variance covariance S 
#                     - matrice d'inertie (correlation) SI
#                     - valeurs propres de S
#                     - vecteurs propres de S
#                     - Les coord des individus et variables  
#                     - Les coord des variables Q-normé dans le cas de l'ACP non normée (le calcul n'est pas bon)
#                     - Inertie total et vecteurs des inerties accumulées.
#                     - La qualité de representation des individus et des variables sur les axes 1 et 2 en %.
#                     - La contribution des individus et des variables sur les axes 1 et 2 en %.
#                     - Le rang de la matrice S: r


ACP = function(df, nf = 2, graphe = FALSE, Scale = TRUE){
  
  ### BLOC DE FONCTIONS ###
  # fonction qui renvoie la moyenne de chaque variables:
  moy_tab = function(x){
    n = nrow(x)
    xbar = (1/n)*t(x)%*%matrix(1,n,1)
    return(xbar)
  }
  #fonction pour l'inverse de la norme d'un vecteur variable:
  norm_var <- function(x) {
    N = diag(t(x)%*%x)^(-1/2)
    return(matrix(N, ncol(x),1))
  }
  #Fonction pour la norme d'un vecteur individu ou un vecteur variable:
  norm_indiv <- function(x) {
    b = rep(NA, nrow(x))
    for (i in seq(1,nrow(x))){
      b[i] = sum(x[i,]^(2))
    }
    return(b)
  }
  #Fonction pour calculer les coord G_k des variables dans le plan(1,2):
  coord = function(lambda, x){
    c = matrix(NA, nrow(x), ncol(x))
    if (length(lambda) == ncol(x)){
      for(i in seq(1:ncol(x))){
        c[,i] = lambda[i]*x[,i]
      }
    }
    return(c)
  }
  #fonction qui permet de centrer notre jeu de donnée 
  centrtab = function(x){         #Les fonction prennent en paramètres des matrices!
    x = as.matrix(x)
    n = nrow(x)
    I = diag(nrow = n)
    A = matrix(1,nrow = n, ncol = n)
    return((I-(1/n)*A)%*%x)
  }
  #fonction qui calcul la matrice de covariance// variance
  #besoin de la métrique sur l'espace des individus D
  covar = function(x){
    x = as.matrix(x)
    n = nrow(x)
    S = (1/n)*t(x)%*%x
    return(S)
  }
  #fonction qui reduit notre matrice 
  reduc = function(x){
    x = as.matrix(x)
    N = diag(diag(covar(x))^(-1/2))
    xcr = x%*%N
    return(xcr)
  }
  
  #on recupère seulement les donées quantitative 
  is_quali = which(!sapply(df, is.numeric))
  if (length(is_quali)>0)
    df = df[,-is_quali]
  df = as.matrix(df)
  #on regarde si dans nos variable il y a des valeurs manquantes NA
  if (any(is.na(df)))
    stop("na values in the dataframe")
  
  #on centre:
  Xc = centrtab(df)
  #Alternative  l'ACP normée et non normée
  if(Scale == TRUE){
    #Cas de l'ACPN, on reduit:
    Xcr = reduc(Xc)}
  else{
    Xcr = Xc #cas de l'acp non normée 
  }
  #on calcule matrice de var/cov
  S = covar(Xcr)
  
  #### Bloc de code pour les métriques ####
  
  nbcol = ncol(df) #nombre de colonnes dans la matrice de variables quantitatives
  if (Scale == TRUE){
    Q = diag(nrow = nbcol)}
  else if(Scale== FALSE){
    Q = diag(diag(S)^(-1)) #matrice d'inverse de variance en diagonale
  }
  
  #on calcule la matrice d'inertie pour le cas de l'ACP non normée :
  SI = S%*%Q 
  
  # On calcul les valeurs/ vecteurs propres:
  v = eigen(S) #renvoi les valeur propres et vecteurs propres
  r = length(v$values) #calcul du rang
  
  #Tracer des valeurs propres:
  barplot(v$values, main = "Inertie des composantes principales",ylab = "valeurs prores")  
  
  #Demande à l'utilisateur d'entrer le nombre d'axes souhaités:
  nf = readline(prompt = "Nombre d'axes: ")
  nf = as.integer(nf)
  if(is.integer(nf) == F){
    stop("la valeur de l'axe doit être numerique.")
  }
  
  #Calcul des facteurs principaux:
  U = Q%*%(v$vectors) #norme inv(Q) 
  
  #on calcul les composantes principales: 
  #Soit la projections des individus dans la nouvelles base
  if (Scale== FALSE){
    Comp_princ = Xcr%*%v$vectors}
  else {
    Comp_princ = Xcr%*%U
  }
  #On recupère seulement les nf axes principaux:
  axes = Comp_princ[, 1:nf]
  
  ##Calcul de tous les G_k:
  var_coord = coord(sqrt(v$values), v$vectors)
 
  #On  les cordonnée des variables Q-othogonal 
  G1 = var_coord[,1]
  G2 = var_coord[,2]
  G3 = var_coord[,3]
 
  #les coordonnées des variables dans le cercle de cor dans le cas de l'acp non normé:
  a1 = matrix(diag(diag(Q^(1/2))%*%t(G1)),nbcol,1)
  a2 = matrix(diag(diag(Q^(1/2))%*%t(G2)), nbcol,1)                     
  a3 = matrix(diag(diag(Q^(1/2))%*%t(G3)), nbcol,1)
  a = cbind(a1,a2,a3)
     ## ALTERNATIVE POUR LE CAS DE LA NON NORMEE 
  
  #### Analyse et interpretation des résultat####
  Ig = sum(v$values) #inertie totale
  I1 = v$values[1]/Ig                 #Calcul de l'inertie sur l'axe 1
  I2 = v$values[2]/Ig                 #Calcul de l'inertie sur l'axe 2
  
  #accumulation des inertie en % 
  accum = cumsum(v$values/sum(v$values)*100) # On arrondi lors de l'appel de la fonction 
  
  
  #### Annalyse indiv ####
  #calcul de la qualité de representation d'un individu sur l'axe1 et l'axe2:  ## REVOIR CTR EN ABS## 
  qlt1=  matrix(diag(norm_indiv(Comp_princ)^(-1)%*%t(axes[,1]^(2))), nrow(df),1)*100 # On exprime en pourcentage
  qlt2 = matrix(diag(norm_indiv(Comp_princ)^(-1)%*%t(axes[,2]^(2))), nrow(df),1)*100
  qlt_indiv= cbind(qlt1,qlt2)
  
  #Calcul de la contribution de l'indiv selon l'axe1 et l'axe2 
  ctr1 = (axes[,1]^(2))/(v$values[1]*nrow(df))*100   # On exprime en pourcentage
  ctr2 = (axes[,2]^(2))/(v$values[2]*nrow(df))*100
  ctr_indiv = cbind(ctr1, ctr2)
  #print(sum(ctr1)) donne 100%
  
  #### Annalyse Variables ####
  #calcul de la qualité de representation d'une variable sur l'axe1 et l'axe2: Contribution relative 
  qlt3 = matrix(diag(norm_indiv(var_coord)^(-1)%*%t(G1^(2))), ncol(df),1)*100
  qlt4 = matrix(diag(norm_indiv(var_coord)^(-1)%*%t(G2^(2))), ncol(df),1)*100
  qlt_var = cbind(qlt3,qlt4)
  
  #Calcul de la contribution de l'indiv selon l'axe1 et l'axe2: 
  ctr3 = (G1^(2))/v$values[1]*100
  ctr4 = (G2^(2))/v$values[2]*100
  ctr_var= cbind(ctr3,ctr4)
  #print(sum(ctr3)) donne 100%
  
  
  #On transforme notre tableau de données centrée reduite en un dataframe
  Xcr = data.frame(Xcr)
  colnames(Xcr) = colnames(df)
  rownames(Xcr) = rownames(df)
  #On transforme la matrice de cov var en dataframe
  S = data.frame(S)
  colnames(S) = colnames(df)
  rownames(S) = colnames(df)
  
  #on convertit en dataframe:
  axes = data.frame(axes)
  rownames(axes) = rownames(df)
  ###
  var_coord = data.frame(var_coord)
  rownames(var_coord) = colnames(df)
  G = var_coord[,1:nf] #on recupère seulement les coord demandées 
  ###
  a = data.frame(a)
  ###
  #On convertit en dataframe nos données:
  qlt_indiv = data.frame(qlt_indiv)
  colnames(qlt_indiv) = c("qlt Axes 1", "qlt Axes 2")
  ###
  ctr_indiv = data.frame(ctr_indiv)
  colnames(ctr_indiv) = c("ctr Axe1", "ctr Axe2")
  ###
  qlt_var = data.frame(qlt_var)
  colnames(qlt_var) = c("qlt Axes 1", "qlt Axes 2")
  ###
  ctr_var =data.frame(ctr_var)
  colnames(ctr_var) = c("ctr Axe1", "ctr Axe2")
  
  ####
  ## On trace les variable dans le cercle de corrélation:
  if (graphe == TRUE & nf == 2){
    layout(matrix(c(1:2),2,1))
    s.label(axes)
    if(Scale == TRUE){
      s.corcircle(G)}
    else{
      s.corcircle(a[,1:2])}
  }
  
  
  res = list(Df = Xcr,mat_varcov = S, mat_inertie = SI, eig = v$values, eig_vec = v$vectors,coor_ind = axes, coor_var = G, coor_var1 = a,Ig = Ig, Inert_cum = accum, qlt_indiv = qlt_indiv, 
             ctr_indiv = ctr_indiv, qlt_var = qlt_var, ctr_var = ctr_var,rank = r )
  
  
}









