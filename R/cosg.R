select_top_n<-function(scores,n_top){
  d <- data.frame(
    x   = data.table::copy(scores),
    indice=seq(1,length(scores)))

  data.table::setDT(d)
  data.table::setorder(d,-x)
  n_top_indice<-d$indice[1:n_top]
  return(n_top_indice)
}

my_cosg<-function(
    object,
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=100,
    group_info = NULL
){

    ### Obtain the cellxgene data
    obj_class = class(object)
    if(obj_class == "Seurat"){
        genexcell<-Seurat::GetAssayData(object = object[[assay]], slot = slot)

        if (groups == 'all'){
            group_info <- Seurat::Idents(object = object)
        }else{
            object <- Seurat::subset(x = object, idents = groups)
            group_info <- Seurat::Idents(object = object)
        }
    }else if(obj_class == "dgCMatrix"){
      genexcell <- object
      group_info = group_info

      }else{
        genexcell <-as(object, "dgCMatrix")
        group_info = group_info

      }





    ### unique groups
    groups_order=sort(unique(group_info))
    n_cluster=length(groups_order)

    if (n_cluster == 1){
        stop('Cannot perform marker gene identification on a single cluster.')}


    n_cell=ncol(genexcell)
    n_gene=nrow(genexcell)
    gene_name=rownames(genexcell)

    ### If sepcifying too many genes to return
    if (n_genes_user>n_gene){
        n_genes_user=n_gene
    }


    cluster_mat=matrix(0,nrow =n_cluster,ncol = n_cell)

    order_i=1
    ### Set gene lambda and gene omega
    for (group_i in groups_order){
        idx_i=group_info==group_i
        cluster_mat[order_i,idx_i]=1
        order_i=order_i+1
    }


    cluster_mat_sparse=as(cluster_mat, "dgCMatrix")
    ### Calculate the cosine similarity
    cosine_sim=proxyC::simil(genexcell,cluster_mat_sparse, method = "cosine",drop0=TRUE)

    pos_nonzero = cosine_sim != 0
    pos_nonzero=which(as.matrix(pos_nonzero),arr.ind = TRUE)

    #### Second-stage
    genexlambda=cosine_sim*cosine_sim
    e_power2_sum=Matrix::rowSums(genexlambda)

    if (mu==1){
         genexlambda[pos_nonzero]=genexlambda[pos_nonzero]/(replicate(ncol(genexlambda),e_power2_sum)[as.matrix(pos_nonzero)])
    }else{
        genexlambda[pos_nonzero]=genexlambda[pos_nonzero]/((
            (1-mu)*genexlambda[pos_nonzero] + mu * (replicate(ncol(genexlambda),e_power2_sum))
        )[as.matrix(pos_nonzero)])
    }

    genexlambda=genexlambda*cosine_sim

    rank_stats_names=data.frame(matrix(matrix(), n_genes_user, length(groups_order),
                        dimnames=list(seq(1,n_genes_user), groups_order)),
                        stringsAsFactors=F)
    rank_stats_scores=data.frame(matrix(matrix(), n_genes_user, length(groups_order),
                        dimnames=list(seq(1,n_genes_user), groups_order)),
                        stringsAsFactors=F)

    order_i=1
    ### Set gene lambda and gene omega
    for (group_i in groups_order){
        idx_i=group_info==group_i
        scores=genexlambda[,order_i]
        global_indices = select_top_n(scores, n_genes_user)
        rank_stats_names[,order_i]=gene_name[global_indices]
        rank_stats_scores[,order_i]=scores[global_indices]

        ### save the group names
        order_i=order_i+1
    }

    colnames(rank_stats_names) <- groups_order
    colnames(rank_stats_scores) <- groups_order

    ###
    ranks_stats=list(
        names=rank_stats_names,
        scores=rank_stats_scores

    )
    ### return
    return(ranks_stats)
}
