############## 
##### HDMD package
##### The 'FactorTransform' function was taken from the HDMD R packge.
##### The packge was archived by CRAN in 2022
##### The code for the function `FactorTransform` was sourced from https://github.com/cran/HDMD
##### All copyrights for this function belog to HDMD package creators under the license GPL (>= 2)

AminoAcids = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
AAbyGroup = c("aliphatic", "cysteine", "acidic", "acidic", "aromatic", "aliphatic", "basic", "aliphatic", "basic", "aliphatic", "aliphatic", "aminic", "proline", "aminic", "basic", "hydroxylated", "hydroxylated", "aliphatic", "aromatic", "aromatic")
AAGroups = c("acidic","aliphatic", "aminic",  "aromatic", "basic","cysteine","hydroxylated", "proline")
small = c(1,2,3,6,12,13,16,17,18)
polar = c(2,3,4,7,9,12,14,15,16,17,19,20)
hydrophobic = c(1, 2,5,6,7,8,9,10,11,17,18,19,20)
#############

# Factor Scores as computed by SAS Factor Analysis of 54 Amino Acid Indices 
AAMetric.Atchley = matrix( 
    c(
        -0.59145974, -1.30209266, -0.7330651,  1.5703918, -0.14550842,
        -1.34267179,  0.46542300, -0.8620345, -1.0200786, -0.25516894,
        1.05015062,  0.30242411, -3.6559147, -0.2590236, -3.24176791,
        1.35733226, -1.45275578,  1.4766610,  0.1129444, -0.83715681,
        -1.00610084, -0.59046634,  1.8909687, -0.3966186,  0.41194139,
        -0.38387987,  1.65201497,  1.3301017,  1.0449765,  2.06385566,
        0.33616543, -0.41662780, -1.6733690, -1.4738898, -0.07772917,
        -1.23936304, -0.54652238,  2.1314349,  0.3931618,  0.81630366,
        1.83146558, -0.56109831,  0.5332237, -0.2771101,  1.64762794,
        -1.01895162, -0.98693471, -1.5046185,  1.2658296, -0.91181195,
        -0.66312569, -1.52353917,  2.2194787, -1.0047207,  1.21181214,
        0.94535614,  0.82846219,  1.2991286, -0.1688162,  0.93339498,
        0.18862522,  2.08084151, -1.6283286,  0.4207004, -1.39177378,
        0.93056541, -0.17926549, -3.0048731, -0.5025910, -1.85303476,
        1.53754853, -0.05472897,  1.5021086,  0.4403185,  2.89744417,
        -0.22788299,  1.39869991, -4.7596375,  0.6701745, -2.64747356,
        -0.03181782,  0.32571153,  2.2134612,  0.9078985,  1.31337035,
        -1.33661279, -0.27854634, -0.5440132,  1.2419935, -1.26225362,
        -0.59533918,  0.00907760,  0.6719274, -2.1275244, -0.18358096,
        0.25999617,  0.82992312,  3.0973596, -0.8380164,  1.51150958
    ), nrow=20, ncol=5, byrow=TRUE)

#Factor Scores as computed by R Factor Analysis (factor.pa.ginv) of 54 Amino Acid Indices   
AAMetric = matrix( 
    c(
        -0.60677211, -1.08940440, -0.92798080,  1.4726202,  0.251417050,
        -1.25436934,  0.46126134, -0.35464046, -1.1040471, -0.115742254,
        1.20123120,  0.24369180, -1.47196724, -0.3243624,  2.127034065,
        1.29519207, -1.46537642, -0.72868426,  0.1482127,  2.078732004,
        -0.95950963, -0.43941703,  1.15581211, -0.5634368, -0.175433772,
        -0.52784144,  1.75873553, -2.14269672,  0.7848519, -0.002451762,
        0.32715696, -0.52087575, -0.42326104, -1.4308910, -0.494428309,
        -1.27633548, -0.59487265,  1.03001883,  0.4256597,  0.014711047,
        1.82248357, -0.53137304, -0.02567446, -0.2432211, -1.571322338,
        -0.97854190, -1.12494822,  0.87095396,  1.4290018, -0.183751797,
        -0.72854024, -1.41768155,  0.30142692, -1.0227830, -0.094615066,
        0.89146155,  1.11359534, -0.85232822, -0.4088754,  0.246465610,
        0.21811094,  1.96849543, -0.15443545,  0.3903353,  0.472209055,
        0.90061235, -0.53619651,  0.02364059, -0.4459318,  0.234557443,
        1.52748331, -0.01654753,  0.65573629,  0.2563726, -2.460536571,
        0.03475833,  0.97046643, -0.84138070,  0.9825769,  0.108360575,
        -0.01042641,  0.65604029,  0.01085433,  0.8716081,  0.093194098,
        -1.32574165, -0.37200927,  1.01019249,  1.4746391,  0.121129607,
        -0.61057537,  0.01937737,  1.50913663, -1.8958892, -0.441521905,
        0.06016331,  0.91703885,  1.35527720, -0.7964405, -0.208006781
    ), nrow=20, ncol=5, byrow=TRUE)

colnames(AAMetric) = c("pah", "pss", "ms", "cc", "ec")
rownames(AAMetric) = AminoAcids

colnames(AAMetric.Atchley) = c("pah", "pss", "ms", "cc", "ec")
rownames(AAMetric.Atchley) = AminoAcids

#####################

Loadings.variation = 
    function(sdev, digits = 5){
        
        var = sdev^2
        
        PTV = var / sum(var)		#proportion var
        CTV = cumsum(PTV)				#cumulative var
        TV = rbind(Lambda=round(var, digits), PTV=round(PTV, digits), CTV=round(CTV, digits))
        TV
    }



FactorTransform = 
    function(Source, Search= AminoAcids, Replace = AAMetric.Atchley, Factor = 1, bycol = TRUE, SeqName = NULL, alignment=FALSE, fillblank=NA){
        
        if(is.matrix(Replace)){
            if(bycol){
                if(ncol(Replace) < Factor) 
                    stop("Column for Factor is greater than size of Replacement Matrix")
                Replace = as.vector(t(Replace[,Factor]))
            }
            else{
                if(nrow(Replace) < Factor) 
                    stop("Row for Factor is greater than size of Replacement Matrix")
                Replace = as.vector(Replace[Factor,])
            }
        }
        
        
        if(is.list(Source)){
            
            if(!missing(SeqName)){
                if(length(SeqName) != length(Source)){
                    print("Length of SeqName not the same as length of Source. Using default.")
                    SeqName = NULL
                }
                else{ SeqNames = SeqName}
            }
            else if (is.null(SeqName)){
                if(is.null(names(Source)))
                    SeqNames = apply(array(1:length(Source)), 1, function(i) paste("Seq",i, sep="")) 			
                else
                    SeqNames = names(Source)		
            }
            
            result = FactorTransform.vector(as.character(Source), Search, Replace)			
            names(result) = SeqNames
        }
        else if(is.vector(Source)){
            if(!missing(SeqName)){
                if(length(SeqName) != length(Source)){
                    print("Length of SeqName not the same as length of Source. Using default.")
                    SeqName = NULL
                }
                else{ SeqNames = SeqName}
            }
            else if (is.null(SeqName)){
                SeqNames = apply(array(1:length(Source)), 1, function(i) paste("Seq",i, sep="")) 
            }
            
            result = FactorTransform.vector(Source, Search, Replace)
            if(length(SeqNames)==1)
                result = list(Seq1=result)
            else
                names(result) = SeqNames
        }
        else if(is.array(Source)){
            result = FactorTransform.vector(Source, Search, Replace)
            rownames(result) = rownames(Source)
            colnames(result) = colnames(Source)	
        }
        else if(is.matrix(Source)){
            result= matrix(nrow= nrow(Source), ncol = ncol(Source))
            for(row in 1:nrow(Source)){
                Seq = Source[row,]
                result[row,] = FactorTransform.default(Seq, Search, Replace)
            }
            rownames(result) = rownames(Source)
            colnames(result) = colnames(Source)
        }
        else{
            stop("Source must be a list, vector, matrix, or array")
        }
        
        if(alignment){	#create matrix
            maxlength = max(sapply(result, length))
            result = t(sapply(result, function(x){append(x, rep(fillblank, times = maxlength-length(x)))}))
            ##		result = matrix(unlist(result), nrow = length(result), byrow = TRUE, dimnames = list(names(result)))
        }
        
        result
    }

FactorTransform.vector = 
    function(Source, Search, Replace){
        
        if(length(Source) ==1){
            Source = unlist(strsplit(Source, NULL, fixed=TRUE))
            SeqMetric = FactorTransform.default(Source, Search, Replace)
        }
        else{
            SeqMetric = lapply(Source, FactorTransform.vector, Search, Replace)			#transform each seq in vector
        }
        SeqMetric
    }

FactorTransform.default=
    function(Source, Search=AminoAcids, Replace=AAMetric.Atchley)
    {
        if(!is.vector(Source))
            stop("Source is not of type vector\n")
        
        if (length(Search) != length(Replace))
            stop("Search and Replace Must Have Equal Number of Items\n")
        
        Changed <- as.character(Source)
        
        for (i in 1:length(Search))
            Changed <- replace(Changed, Changed == Search[i], Replace[i])
        
        as.numeric(Changed)
    }
## adapted from Marc Schwartz code at https://stat.ethz.ch/pipermail/r-help/2006-July/108829.html
