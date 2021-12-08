############ Normalization ###################
# This function is used for performing normalization on gene expression data.
Normalization <- function(matrix, option,meta_ref){
    # CPM TPM QN
    if ( option == "CPM"){

        if(class(matrix)[[1]] == "dgCMatrix"){
            matrix <- t(t(matrix)/colSums2(matrix))*1e6
        }else{
            matrix <- cpm(matrix)
        }

    }else if (option=="TPM"){
        if(length(grep("ENSG000",rownames(matrix))) > 100){
            data(db)
            # bulkdata2 <- bulkdata[rowSums(bulkdata > 1) >= 230, ]  ##  detectable if at least two cells contain more than 1 transcript from the gene
            geneNames <- db[match(rownames(matrix), db[,"ensembl_gene_id"]), ]
            geneNames2 <- geneNames[!(is.na(geneNames[,"ensembl_gene_id"] != "")), ]
            geneNames2$length <- geneNames2$end_position - geneNames2$start_position
            matrix <- matrix[geneNames2$ensembl_gene_id, ]


            matrix <- tpmnorm(matrix, geneNames2$length )
        }else{
            data(human_lengths)
            rownames(matrix) <- tolower(rownames(matrix))
            names(human_lengths)  <- tolower(names(human_lengths))

            A  <- intersect(rownames(matrix),names(human_lengths))
            matrix  <- matrix[A,]
            human_lengths  <- human_lengths[A]
            rate  <- matrix / human_lengths
            apply(rate,2,function(x) 1e6*x/sum(x))
            rownames(matrix)  <- toupper(rownames(matrix))
        }

        ####################################################################################


    }else if(option=="QN"){

        matrix_rownames <- rownames(matrix)
        matrix_colnames <- colnames(matrix)

        matrix = normalize.quantiles(as.matrix(matrix))

        rownames(matrix) <- matrix_rownames; colnames(matrix) <- matrix_colnames
    }else if (option=="TMM"){# CPM counts coming from TMM-normalized library sizes; https://support.bioconductor.org/p/114798/

        matrix <- calcNormFactors(matrix, method = "TMM")
        if(class(matrix)[[1]] == "dgCMatrix"){
            matrix <- t(t(matrix)/colSums2(matrix))*1e6
        }else{
            matrix <- cpm(matrix, log=FALSE)
        }

    }

    #Seurat pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000), Scran (centre_size_factors=TRUE) and SCnorm (K=1, conditions=rep(c(1), each=1000))

    return(matrix)

}

########### get_input_ensemble ###########
get_input_ensemble <- function(count_bulk, ref_matrix, meta_bulk, meta_ref, true_frac = NULL,params) {

    data_type <- params$data_type

    # Subset common genes
    gene <- intersect(rownames(ref_matrix), rownames(count_bulk))

    count_bulk <- count_bulk[gene,]
    ref_matrix <- ref_matrix[pmatch(gene,rownames(ref_matrix)),]
    tmp_count_bulk <- count_bulk

    # Normalization
    if(params$TNormalization == "TMM"){
        ct = cbind(ref_matrix,count_bulk)
        ct <- DGEList(counts = ct)
        ct = Normalization(ct,"TMM")
        ref_matrix = ct[,1:ncol(ref_matrix)]
        count_bulk = ct[,(ncol(ref_matrix)+1):ncol(ct)]
        rm(ct)
    }else if(params$TNormalization != "none"){
        count_bulk <- Normalization(count_bulk,params$TNormalization)
    }else if(params$CNormalization != "none"){
        ref_matrix <- Normalization(ref_matrix,params$CNormalization)
    }



    if(params$data_type %in% c("singlecell-rna","rna-seq")){
        tmp_count_bulk <- Normalization(tmp_count_bulk,"TPM")
        count_bulk <- count_bulk[rownames(tmp_count_bulk),]
        gene <- intersect(rownames(ref_matrix), rownames(count_bulk))

        count_bulk <- count_bulk[gene,]
        ref_matrix <- ref_matrix[pmatch(gene,rownames(ref_matrix)),]
        tmp_count_bulk <- tmp_count_bulk[gene,]
        }

    # Scaling
    if(params$Scale == "log"){
        count_bulk = log2(count_bulk+1)
        if(class(ref_matrix)[[1]] == "dgCMatrix"){
            ref_matrix@x <- log2(ref_matrix@x + 1)
        }else{
            ref_matrix <- log2(ref_matrix+1)
        }
    }

    data_set <- get_set_basis(count_bulk = count_bulk,ref_matrix = ref_matrix,meta_bulk = meta_bulk,meta_ref = meta_ref,params = params,true_frac = true_frac)
    data_set <- list(data_set)
    names(data_set) <- params$data_name


    return(data_set)
}


###################### get_set_basis #################
get_set_basis <- function(count_bulk = count_bulk,ref_matrix = ref_matrix,meta_bulk = meta_bulk,meta_ref = meta_ref,params = params,true_frac = true_frac){

    data_set <- list()
    data_type <- params$data_type

    if(data_type == "singlecell-rna"){
        data_set$data$data_c <- ref_matrix
        data_set$data$meta_ref <- meta_ref

        data_set$data$data_t <- count_bulk

        cell_type <- sort(unique(meta_ref$deconv_clust))
        K <- length(cell_type)
        data_set$annotation$pure <- t(sapply(meta_ref$deconv_clust, function(x) cell_type == x))
        rownames(data_set$annotation$pure) <- rownames(meta_ref)
        colnames(data_set$annotation$pure) <- cell_type
        data_set$annotation$pure_samples <- lapply(1:K, function(x){which(data_set$annotation$pure[,x])})
        names(data_set$annotation$pure_samples) <- cell_type
        data_set$annotation$prep <- FALSE
        if(is.null(true_frac)){
            true_frac <- matrix(NA, nrow = ncol(count_bulk), ncol = K)
        }  else{
            true_frac <- true_frac[,cell_type]
            true_frac <- sum_to_one(true_frac)
        }
        data_set$annotation$mixture <- rbind(1*(data_set$annotation$pure), true_frac)

        data_set$annotation$data_type <- data_type
        data_set$name <- params$data_name

        data_set$notes <- params

        data_set$ensemble$test_bulk <- count_bulk
    }else{
        data_set$data$data_t <- count_bulk
        data_set$data$meta_ref <- meta_ref
        data_set$data$data_c <-ref_matrix


        if(is.null(meta_ref)){
            cell_type <- colnames(ref_matrix)

            K <- length(cell_type)
            data_set$annotation$pure <- t(sapply(colnames(ref_matrix), function(x) cell_type == x))
            rownames(data_set$annotation$pure) <- colnames(ref_matrix)
        }else{
            cell_type <- sort(unique(meta_ref$deconv_clust))
            K <- length(cell_type)
            data_set$annotation$pure <- t(sapply(meta_ref$deconv_clust, function(x) cell_type == x))
            rownames(data_set$annotation$pure) <- rownames(meta_ref)
        }


        colnames(data_set$annotation$pure) <- cell_type

        data_set$annotation$pure_samples <- lapply(1:K, function(x){which(data_set$annotation$pure[,x])})
        names(data_set$annotation$pure_samples) <- cell_type

        data_set$annotation$prep <-FALSE

        if(is.null(true_frac)){
            true_frac <- matrix(NA, nrow = ncol(count_bulk), ncol = K)
        }  else{
            true_frac <- true_frac[,cell_type]
            true_frac <- sum_to_one(true_frac)
        }

        data_set$annotation$mixture <- rbind(1*(data_set$annotation$pure), true_frac)
        data_set$annotation$data_type <- data_type
        data_set$name <- params$data_name
        data_set$notes <- params
        data_set$ensemble$test_bulk <- count_bulk
    }

    return(data_set)
}

############## cv_set ############
# Create cross validation dataset
# cv_set <- function(x,count_bulk = count_bulk,ref_matrix = ref_matrix,meta_bulk = meta_bulk,meta_ref = meta_ref,params = params,true_frac = true_frac){
# 
#     data_type <- params$data_type
# 
#     picked <-rownames(count_bulk)[x]
#     test_bulk <-count_bulk[pmatch(picked,rownames(count_bulk)),]
#     train_bulk <-count_bulk[-pmatch(picked,rownames(count_bulk)),]
# 
#     train_ref <- ref_matrix[-pmatch(picked,rownames(ref_matrix)),]
# 
#     data_set <- list()
# 
#     if(data_type == "singlecell-rna"){
#         data_set$data$data_c <- train_ref
#         data_set$data$meta_ref <- meta_ref
#         data_set$data$data_t <- train_bulk
# 
#         cell_type <- sort(unique(meta_ref$deconv_clust))
#         K <- length(cell_type)
#         data_set$annotation$pure <- t(sapply(meta_ref$deconv_clust, function(x) cell_type == x))
#         rownames(data_set$annotation$pure) <- rownames(meta_ref)
#         colnames(data_set$annotation$pure) <- cell_type
#         data_set$annotation$pure_samples <- lapply(1:K, function(x){which(data_set$annotation$pure[,x])})
#         names(data_set$annotation$pure_samples) <- cell_type
#         data_set$annotation$prep <-FALSE
#         if(is.null(true_frac)) true_frac <- matrix(NA, nrow = nrow(meta_bulk), ncol = K) else true_frac = true_frac[,cell_type]
#         data_set$annotation$mixture <- rbind(1*(data_set$annotation$pure), true_frac)
# 
#         data_set$annotation$data_type <- params$data_type
#         data_set$name <- params$data_name
# 
#         data_set$notes <- params
# 
#         data_set$ensemble$test_bulk <- test_bulk
# 
# 
#     }else{
# 
#         data_set$data$data_t <- train_bulk
# 
#         data_set$data$data_c <-train_ref
# 
#         data_set$data$meta_ref <- meta_ref
# 
# 
#         if(is.null(meta_ref)){
#             cell_type <- colnames(ref_matrix)
#             K <- length(cell_type)
#             data_set$annotation$pure <- t(sapply(colnames(ref_matrix), function(x) cell_type == x))
#             rownames(data_set$annotation$pure) <- colnames(ref_matrix)
#         }else{
#             cell_type <- sort(unique(meta_ref$deconv_clust))
#             K <- length(cell_type)
#             data_set$annotation$pure <- t(sapply(meta_ref$deconv_clust, function(x) cell_type == x))
#             rownames(data_set$annotation$pure) <- rownames(meta_ref)
#         }
# 
#         colnames(data_set$annotation$pure) <- cell_type
# 
#         data_set$annotation$pure_samples <- lapply(1:K, function(x){which(data_set$annotation$pure[,x])})
#         names(data_set$annotation$pure_samples) <- cell_type
# 
#         data_set$annotation$prep <-FALSE
# 
#         if(is.null(true_frac)) true_frac <- matrix(NA, nrow = ncol(count_bulk), ncol = K) else true_frac = true_frac[,cell_type]
# 
#         data_set$annotation$mixture <- rbind(1*(data_set$annotation$pure), true_frac)
# 
#         data_set$annotation$data_type <- params$data_type
# 
#         data_set$name <- params$data_name
# 
# 
#         data_set$notes <- params
# 
#         data_set$ensemble$test_bulk <- test_bulk
# 
#     }
# 
#     return(data_set)
# }

######## tpmnorm ###########
tpmnorm <- function(counts,len) {
    x <- counts/len
    return(t(t(x)*1e6/colSums2(x)))
}

########## rm_zero ###########
#remove matrix column that contains all 0s

rm_zero = function(mat,type = "Matrix"){

    if(type == "dgCMatrix"){
        mat <- mat[rowSums2(mat)>0,]
        mat <- mat[,colSums2(mat)>0]
        # mat = mat[!rowVars(mat) == 0,]
    }else{
        mat <- mat[rowSums(mat)>0,]
        mat <- mat[,colSums(mat)>0]
        # mat = mat[!apply(mat, 1, function(x) var(x) == 0),]
    }

    return(mat)
}

############### sum_to_one ###############
sum_to_one <- function(matrix){

    matrix[matrix  < 0] <- 0
    matrix <- apply(matrix,2,function(x){
        x[is.na(x)] = 0
        return(x)})

    if(any(rowSums(matrix) == 0)){
        matrix[rowSums(matrix) == 0,] <- 1/ncol(matrix)
    }
    matrix <- matrix/rowSums(matrix)
}


ord_name <- function(df,df_true){

    df <- df[,colnames(df_true)]
    df <- df[rownames(df_true),]
    return(df)
}
############## filterzerovar ##############3
filterzerovar <- function(mat){

    if(class(mat)[[1]] == "dgCMatrix"){
        mat <- mat[!rowVars(mat) == 0,]
    }else{
        mat <- mat[!apply(mat, 1, function(x) var(x) == 0),]
    }
    return(mat)
}


########### Batch Correction #########
########### B mode ###########
B_mode <- function(mix, ref,scale,phat){
    if(scale == "linear"){
        gene <- intersect(rownames(mix),rownames(ref))
        mix <- mix[gene,]
        ref <- ref[gene,]
        K <- ncol(ref)

        mstar <- ref %*% t(phat)
        lmstar <- log2(mstar+1)
        lmix <- log2(mix+1)

        data <- cbind(lmix,lmstar)
        batch <- c(rep(1,ncol(lmstar)),rep(2,ncol(lmstar)))
        combat_edata1 <- ComBat(dat=data, batch=batch)

        cmix <- 2^(combat_edata1[,1:ncol(lmix)])-1
    }else{
        gene <- intersect(rownames(mix),rownames(ref))
        mix <- mix[gene,]
        ref <- ref[gene,]
        K <- ncol(ref)
        mstar <- ref %*% t(phat)

        data <- cbind(mix,mstar)
        batch <- c(rep(1,ncol(mstar)),rep(2,ncol(mstar)))
        combat_edata1 <- ComBat(dat=data, batch=batch)

        cmix <- combat_edata1[,1:ncol(mix)]
    }


    return(cmix)
}


########### S mode ###########
S_mode <- function(sig_matrix,mix,meta_ref,ref_matrix,transformation){
    ref_matrix <- as.matrix(ref_matrix[rownames(sig_matrix),])
    ct <- unique(meta_ref$deconv_clust)
    c <- length(ct)
    mu <- table(meta_ref$deconv_clust)/length(meta_ref$deconv_clust)
    #F_mtx = matrix(0,nrow = ncol(mix),ncol = c)
    set.seed(2021)
    F_star <- sapply(1:c, function(x) rnorm(ncol(mix),mean = mu[x],sd = 2*mu[x]))
    F_star[F_star<0] =0
    F_star <- F_star/rowSums(F_star)

    if(transformation == "TPM"){
        tmp <- matrix(0,nrow = nrow(ref_matrix),ncol = c)
        res <- matrix(0,nrow = nrow(ref_matrix),ncol = ncol(mix))
        for (j in 1:ncol(mix)) {
            for (i in 1:c) {
                tmp[,i] = colSums(sample_frac(as.data.frame(t(ref_matrix[,which(meta_ref$deconv_clust == ct[i])])),F_star[j,i]))
            }
            res[,j] = rowSums(tmp)

        }

    }

    data = cbind(lmix,lmstar)
    batch = c(rep(1,ncol(lmstar)),rep(2,ncol(lmstar)))
    combat_edata1 = ComBat(dat=data, batch=batch)

}



get_os <- function(){
    os_type =.Platform$OS.type
    if(os_type == "unix"){
        sysinf <- Sys.info()
        if (!is.null(sysinf)){
            os <- sysinf['sysname']
            if (os == 'Darwin')
                os <- "osx"
        } else { ## mystery machine
            os <- .Platform$OS.type
            if (grepl("^darwin", R.version$os))
                os <- "osx"
            if (grepl("linux-gnu", R.version$os))
                os <- "linux"
        }
        tolower(os)
    }else{
        os = os_type
    }
    return(os)
    
}
