kms_cv <- function(data,k = 5,centers = NULL,max_center = 20,alpha = .1,seed = NULL){
    
    if(!is.null(seed)){
        set.seed(seed)
    }
    
    data <- parsnip::maybe_data_frame(data)
    
    column_order <- base::colnames(data)
    
    numeric_col <- base::sapply(data,base::is.numeric)
    all_col <- base::ncol(data)
    
    if(sum(numeric_col) != sum(all_col)){
        
        base::warning(base::sprintf('Dropping %s columns, numerics are used only.',sum(all_col) - sum(numeric_col)))
        
        data_ <- data[,numeric_col]
        
    }else{
        data_ <- data
    }
    
    if(is.null(centers)){
        
        baseline_tot_withinss <- stats::kmeans(data_,centers = 1)[['tot.withinss']]
        
        scaled_inertias <- c()
        
        for(i in 2:max_center){
            
            trial <- stats::kmeans(data_,centers = i)
            
            scaled_inertias[i] <- (trial[['tot.withinss']] / baseline_tot_withinss) + alpha*i
            
        }
        
        centers <- base::which.min(scaled_inertias)
    }
    
    clustered <- stats::kmeans(data_,centers = centers)[['cluster']]
    
    data_[['cls']] <- clustered
    
    current_prop <- base::prop.table(base::table(data_[['cls']]))
    
    current_prop_df <- dplyr::rename(parsnip::maybe_data_frame(current_prop),'cls' = 'Var1','Prop' = 'Freq')
    
    current_prop_df[['cls']] <- base::as.integer(current_prop_df[['cls']])
    
    ideal_sample_size <- base::round(base::nrow(data_) / k)
    
    ## Resampling
    
    resamples <- list()
    
    for(i in 1:k){
        
        
       props <- current_prop_df[['Prop']]
        
       group_size <- ideal_sample_size * props
        
       resample <- purrr::map2_dfr(dplyr::group_split(data_,cls),group_size, ~ slice_sample(.x, n = .y))
        
       resample[['id']] <- base::sprintf('Fold%s',i)
        
       resamples[[i]] <- resample
        
       data_ <- dplyr::anti_join(data_,resample,by = base::colnames(data_)) 
        
    }
    
    #Slacks
    for(i in 1:base::nrow(data_)){
        
        random_fold <- base::sample(1:k,size = 1)
        
        to_kick <- data_[1,]
        
        to_kick[['id']] <- base::sprintf('Fold%s',random_fold)
        
        resamples[[i]] <- parsnip::maybe_data_frame(base::rbind(resamples[[i]],to_kick))
        
        data_ <- dplyr::anti_join(data_,to_kick,by = base::colnames(data_))
        
    }
    
    resample_union <- base::do.call(base::rbind,resamples)
    
    resample_union[['cls']] <- NULL
    
    join_cols <- base::setdiff(base::colnames(resample_union),'id')
    
    to_nest <- dplyr::left_join(resample_union,data,by=join_cols)
    
    order_cols <- c(column_order,'id')
    
    to_nest <- dplyr::select(to_nest,dplyr::all_of(order_cols))
    
    tidyr::nest(dplyr::group_by(to_nest,id),.key = 'splits')
}