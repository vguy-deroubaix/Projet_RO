
#= NOM1 - Prénom 1
   NOM2 - Prénom 2
   N'oubliez pas de modifier ce commentaire, ainsi que le nom du fichier! =#

   using JuMP, GLPK


   # Implementation 1 Miller-Tucker-Zemlin
   function resolutionRapideDuTSP(C::Matrix{Int64})
   
        m::Model = Model(GLPK.Optimizer)
   
        nbSite::Int64 = size(C, 1)
   
        vectorSite::Vector{Int64} = collect(1:nbSite)
   
        @variable(m, nbSite >= t[2:nbSite] >=0 , Int)
        @variable(m, x[i = 1:nbSite, j = 1:nbSite], Bin)
    
        @objective(m, Min, sum(sum(x[i,j]*C[i,j] for i in 1:nbSite) for j in 1:nbSite))
   
        @constraint(m, contrA[i in 1:nbSite], sum(x[i,j] for j in deleteat!(copy(vectorSite), i)) == 1)
   
        @constraint(m, contrB[j in 1:nbSite], sum(x[i,j] for i in deleteat!(copy(vectorSite), j)) == 1)
    
        @constraint(m, contrC[i in 2:nbSite, j in 2:nbSite], t[i]-t[j] + nbSite*x[i,j] <= nbSite - 1) 
   
        optimize!(m)
   
        # Affichage des résultats
        status = termination_status(m)
   
        if status == MOI.OPTIMAL
           println("Problème résolu à l'optimalité")
           println("z = ",objective_value(m)) # affichage de la valeur optimale
           println("x = ",value.(x))
        elseif status == MOI.INFEASIBLE
           println("Problème impossible!")
        elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
           println("Problème non borné!")
        end
        return m
   end
   
   
   # Implementation 2 Dantzig-Fulkerson-Johnson
   
   # Modèle de problème d'affectation 
    function resolutionPbAffectation(C::Matrix{Int64})
   
        m::Model = Model(GLPK.Optimizer)
   
        nbSite::Int64 = size(C,1)
        vectorSite::Vector{Int64} = collect(1:nbSite)
   
        @variable(m, x[i = 1:nbSite, j = 1:nbSite], Bin)
        @objective(m, Min, sum(sum(x[i,j]*C[i,j] for i in 1:nbSite) for j in 1:nbSite))
   
        # contraintes de 1 et 2
        @constraint(m, contrA[i in 1:nbSite], sum(x[i,j] for j in deleteat!(copy(vectorSite), i)) == 1)
        @constraint(m, contrB[j in 1:nbSite], sum(x[i,j] for i in deleteat!(copy(vectorSite), j)) == 1)
   
        optimize!(m)
   
        return m
   end
   
    function resolutionMalineDuTSP(C::Matrix{Int64})
   
        model = resolutionPbAffectation(C)

        X::Matrix{Float64} = value.(model[:x])
       
        cycleMin::Vector{Tuple} = minCycle(X)# recupérer les résultats => deduire le plus petit sous tour 
       
        while size(cycleMin,1) != size(X,1)# Tant que l'on a des sous tours :
       
            @constraint(model,sum(model([:x])[i,j] for (i,j) in cycleMin) <= 1)# resoudre le pb d'affectation avec contrainte supplémentaire : casser le plus petit sous tour  ( xij + xji <= 1)   addconstraint ?
            model = resolutionMalineDuTSP(C)
            X = value.(model[:x]) 
            cycleMin = minCycle(X)

        end# fin tant que 
   
        # afficher résultats 
        status = termination_status(model)
   
        if status == MOI.OPTIMAL
           println("Problème résolu à l'optimalité")
           println("z = ",objective_value(model)) # affichage de la valeur optimale
           println("x = ",value.(model[:x]))
        elseif status == MOI.INFEASIBLE
           println("Problème impossible!")
        elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
           println("Problème non borné!")
        end
           
    end

    function minCycle(X::Matrix{Float64})

        A::Matrix{Tuple{Int64,Int64}} = [] # Matrice dans laquelle on stock les cycles
        min::Vector{Tuple{Int64,Int64}} 
        
        #vecteur des indices i et j que l'on racoursira au fur et a mesure de la boucle 
        height::Vector{Int64} = collect(1:size(X,1)) 
        witdh::Vector{Int64}  = collect(1:size(X,2))
        
        #indice i et j qui sert a voyager dans X
        i::Int64 = height[1]
        j::Int64 = width[1] 
        k::Int64 = 1 # k sert d'indice dans A
        l::Int64 = 1 # l sert a traverser le tableau d'indice de j
        cyclesuivant::Bool = false

        #le but est de parcurir le moins de d'element possible normalement on est sur du n(n+1)/2
        while !isempty(height)
            while !isempty(width)
                if X[i,j] == 1 
                    #=
                    #tant que on a pas fini de mettre dans le table le cycle 
                    qu'on est en train de regarder on change pas de cycle
                    =#
                    if !cyclesuivant 
                        push!(A[k],(i,j))
                    else 
                        k = k+1
                        push!(A[k,(i,j)])
                    end 
                    #=
                    On supprime les indices qu'on a déjà regarde car il ne peut y avoir qu'une seule
                    variable a 1 sur ça ligne et ça colone 
                    donc si on trouve une variable a 1 on a plus besion de regarde ni ça ligne ni ça colone 
                    =#
                    delete!(height,i)
                    delete!(width,j)
                    #=
                    A la fin on fait un tour de boucle alors que les tableau d'indice sont vide 
                    donc faut verifier que les tableau soit pas vide pour pas fait un segfault
                    =#
                    if !isempty(height) && !isempty(witdh)
                        #=
                        si j est dans witdh alors on a pas fini le cycle 
                        =#
                        if in(j,height) 
                            i = j
                        else 
                            #=
                            si on a fini le cycle on revient au début du tableau d'indice
                            pour pouvoir trouver le début du prochain cycle
                            =#
                            i = height[1]
                            cyclesuivant = true
                        end 
                        #=
                        on revient ensuite au début des j pour pouvoir traverser toute la ligne 
                        en sautant les colones que l'on a déjà visiter
                        =#
                        j = width[1]
                        l = 1
                    end
                else 
                    j = width[l+1]
                end
            end
        end
        min = A[1]
        for i in 1:size(A,1)
            if size(A[i],1) < size(min,1)
                min = A[i]
            end
        end
        return min
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