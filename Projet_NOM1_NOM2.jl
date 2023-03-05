
#= NOM1 - Prénom 1
   NOM2 - Prénom 2
   N'oubliez pas de modifier ce commentaire, ainsi que le nom du fichier! =#

using JuMP, GLPK

# Une fonction qu'on devra écrire (en fait, il y en aura une pour chaque méthode, plus d'autres fonctions utiles...)
function resolutionRapideDuTSP(C::Matrix{Int64})

    #TODO
    m::Model = Model(GLPK.Optimizer)

    nbSite::Int64 = size(C, 1)

    vectorSite::Vector{Int64} = collect(1:nbSite-1)

    @variable(m, nbSite >= t[1:nbSite] >=0 , Int)
    @variable(m, x[i = 1:nbSite, j = 1:nbSite], Bin)

    @objective(m, Min, sum(sum(x[i,j]*C[i,j] for i in 1:nbSite) for j in 1:nbSite))

    @constraint(m, contrA[i in 1:nbSite], sum(x[i,j] for j in deleteat!(vectorSite, i)) == 1)

    @constraint(m, contrB[j in 1:nbSite], sum(x[i,j] for i in deleteat!(vectorSite, j)) == 1)
 
    @constraint(m, contrC[i in 1:nbSite, j in 1:nbSite,], t[i]-t[j] + nbSite*x[i,j] <= nbSite - 1) 
end



# fonction qui prend en paramètre un fichier contenant un distancier et qui retourne le tableau bidimensionnel correspondant
function parseTSP(nomFichier::String)
    # Ouverture d'un fichier en lecture
    f::IOStream = open(nomFichier,"r")

    # Lecture de la première ligne pour connaître la taille n du problème
    s::String = readline(f) # lecture d'une ligne et stockage dans une chaîne de caractères
    tab::Vector{Int64} = parse.(Int64,split(s," ",keepempty = false)) # Segmentation de la ligne en plusieurs entiers, à stocker dans un tableau (qui ne contient ici qu'un entier)
    n::Int64 = tab[1]

    # Allocation mémoire pour le distancier
    C = Matrix{Int64}(undef,n,n)

    # Lecture du distancier
    for i in 1:n
        s = readline(f)
        tab = parse.(Int64,split(s," ",keepempty = false))
        for j in 1:n
            C[i,j] = tab[j]
        end
    end

    # Fermeture du fichier
    close(f)

    # Retour de la matrice de coûts
    return C
end




#= Exemple de script, qui résout ici les instances jusqu'à une taille de 40. Ce script devra être modifié/adapté suivant les besoins. 
    La macro @time mesure le temps (elapsed) d'exécution d'une fonction
    La consommation mémoire est aussi indiquée mais la valeur indiquée n'est correcte que pour un code 100% codé en Julia (et GLPK est codé en C) =#
function ExempleScriptTSP()
    # Première exécution sur l'exemple pour forcer la compilation si elle n'a pas encore été exécutée
    C::Matrix{Int64} = parseTSP("plat/exemple.dat")

    resolutionRapideDuTSP(C)


    file::String = ""
    # Série d'exécution avec mesure du temps pour des instances symétriques
    for i in 10:10:40
        file = "plat/plat$i.dat"
        C = parseTSP(file)
        println("Instance à résoudre : plat$i.dat")
        @time resolutionRapideDuTSP(C)
    end

    # Série d'exécution avec mesure du temps pour des instances asymétriques
    for i in 10:10:40
        file = "relief/relief$i.dat"
        println("Instance à résoudre : relief$i.dat")
        C = parseTSP(file)
        @time resolutionRapideDuTSP(C)
    end
end