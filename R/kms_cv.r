kms_cv <- function(data,k = 5,centers = NULL,max_center = 20,alpha = .1,seed = NULL){
    
    if(!is.null(seed)){
        set.seed(seed)
    }
    
    data <- parsnip::maybe_data_frame(data)
    
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
    
    
}