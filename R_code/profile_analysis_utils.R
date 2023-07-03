library(dplyr)

# profiles_dict are pairs of name  and filename
# the profiles in filenames will be loaded into new variables by the names(files_dict)
read_profiles_and_assign_vars <- function (profiles_dict, dir_name){
    keys<- names(files_dict)
    files<- files_dict
    for (i in 1:length(keys)) {
        full_path_profile <- file.path(profiles_dir, files[i])
        profile <- read.csv(full_path_profile, na.strings = "", row.names = 1)
        assign(keys[i],   profile, envir = .GlobalEnv)
    }
}


invisible(capture.output({
    library(matrixStats)
    library(dplyr)
}))
cutoff_complete_cases <- function(cutoff,profile_df) {
    profile_df[abs(profile_df) < cutoff] <- NA
    data.frame(profile_df[complete.cases(profile_df), ])
}

cutoff_scale_complete_cases <- function(cutoff, profile_df) {
    df <- data.frame()
    df <- cutoff_complete_cases(cutoff, profile_df) %>%
        scale() %>% as.data.frame() %>% mutate_all(~as.numeric(.))
    df
}

ALL_ROWS = -1
get_top_sds <- function(sds_count, df) {
    # print(sds_count)
    # print(nrow(df))
    sds <- rowSds(as.matrix(df), na.rm = TRUE)
    o <- order(sds, decreasing = TRUE)
    if (sds_count > nrow(df)) {
        return (NA)
    }
    if (sds_count == ALL_ROWS) {
        return(df[o,])
    }
    df[o[1:sds_count],]
}



cutoff_and_get_top_sds <- function(sds_count, cutoff, profile_df) {
    prf_cutoffed_df <- cutoff_scale_complete_cases(cutoff, profile_df)
    return (get_top_sds(sds_count, prf_cutoffed_df))
}


get_function_grid_old <- function(top_sds, cutoffs, df, my_func, method = NULL) {
    # print(method)
    # result <- outer(top_sds, cutoffs, Vectorize(function(a, b) my_func(a, b, df)))
    result_mat <- matrix(nrow = 0, ncol = 0)
    if (is.null(method)) {
        result_mat <- outer(top_sds, cutoffs, Vectorize(function(a, b) my_func(a, b, df)))
    }else{
        result_mat <- outer(top_sds, cutoffs, Vectorize(function(a, b) my_func(a, b, df, method)))
    }

    result_df <- data.frame(result_mat)
    top_sds[length(top_sds)] = "ALL_ROWS"
    rownames(result_df) <- top_sds
    colnames(result_df) <- cutoffs
    result_df
}


cutoff_and_top_sds_list_exec_func <- function(sds_count_list, cutoff, profile_df, my_func, method="complete"){
    prf_cutoffed_df <- cutoff_scale_complete_cases(cutoff, profile_df)
    cleaned_scaled_prf_lst <- get_top_sds_list(sds_count_list, prf_cutoffed_df)
    results = list()
    for (mat in cleaned_scaled_prf_lst){
        if (any(is.na(mat))) {
            results <- append(results,NA)
        } else {
            pci_result = my_func(mat, method)
            results <- append(results, pci_result)
        }
    }
    results
}

get_function_grid <- function(top_sds, cutoffs, m_df, m_func, m_method = NULL) {
    result_mat <- matrix(nrow = 0, ncol = 0)
    if (is.null(m_method)) {
        results_in_lists <- sapply(cutoffs,
                function(x, std_list=top_sds, df= m_df, my_func=m_func){
                    cutoff_and_top_sds_list_exec_func(std_list, x, df, my_func, "complete")
        })
        result_mat = rbind(results_in_lists)
    } else {
        results_in_lists <- sapply(cutoffs,
                function(x, std_list=top_sds, df= m_df, my_func=m_func, method=m_method){
                    cutoff_and_top_sds_list_exec_func(std_list, x, df, my_func, method)
        })
        result_mat = rbind(results_in_lists)
    }

    result_df <- data.frame(result_mat)
    top_sds[length(top_sds)] = "ALL_ROWS"
    rownames(result_df) <- top_sds
    colnames(result_df) <- cutoffs
    result_df
}


######################## FOR CLASSIFICATION EVAL ######################################33
get_top_sds_list <- function(sds_count_list, df) {
    sds <- rowSds(as.matrix(df), na.rm = TRUE)
    o <- order(sds, decreasing = TRUE)
    results = list()
    for (sds_count in sds_count_list) {
        if (sds_count == ALL_ROWS){
            sds_count = length(o)
        }
        if (sds_count > nrow(df)){
            results <- append(results, NA)
        } else {
            sds_mat <- df[o[1:sds_count],]
            results <- append(results, list(sds_mat))
        }
    }
    return (results)
}



avg_func_res_per_df_new <- function(dfs_list, m_top_sds, names,  m_func, m_method, m_cutoff){
    result_vector = sapply(dfs_list, function (x, method= m_method, top_sds=m_top_sds, cutoff=m_cutoff, my_func=m_func) {
        pci_list <-  cutoff_and_top_sds_list_exec_func(top_sds, cutoff, x, my_func, method)
        # pci_list <- cutoff_and_top_sds_list_get_PCI(top_sds, cutoff, x, method)
        as.numeric(mean(as.numeric(pci_list), na.rm = TRUE) )
    })
    names(result_vector) <- names
    return (result_vector)
}

avg_func_res_sds_df_new <- function(dfs_list, top_sds, names, my_func, cutoff){
    avg_pci_sds             <- avg_func_res_per_df_new(dfs_list, top_sds, names, my_func, "complete", cutoff)
    avg_pci_sds_avg         <- avg_func_res_per_df_new(dfs_list, top_sds, names, my_func, "average",  cutoff)
    avg_pci_sds_ward        <- avg_func_res_per_df_new(dfs_list, top_sds, names, my_func, "ward.D2",   cutoff)

    pci_sds_df  <- data.frame(
        "complete"  = avg_pci_sds,
        "average"   = avg_pci_sds_avg,
        "ward"      = avg_pci_sds_ward)
    as.data.frame(round(t(pci_sds_df),3))
}



# --------------------------------  calculation on a grid +  N  -------------------------
avg_func_n_vals <- function(dfs_list, m_N, names, my_func){
    result_vector = sapply(dfs_list, function (x, N=m_N) {
        my_func(x, N)
    })
    names(result_vector)<- names
    return (result_vector)
}

avg_func_n_df <- function(dfs_list, N, names, title, my_func){
    avg_vals_max_n    <- avg_func_n_vals(dfs_list, N, names, my_func)
    vals_sds_df  <- data.frame(avg_vals_max_n)
    colnames(vals_sds_df) <- title
    as.data.frame(round(t(vals_sds_df),3))
}


#------------------------------- Classification Eval - Average of max N in grid

extract_avg_of_n_max_df <- function(df,N) {
    to_vector <- as.numeric(unlist(df))
    sorted_vector <- sort(to_vector, decreasing = TRUE)
    if (N > length(sorted_vector)){N = length(sorted_vector)}
    na.rm = TRUE
    sorted_values <- as.numeric(sorted_vector[1:N], na.rm= TRUE)
    mean(sorted_values)
}


avg_vals_max_n_df <- function(dfs_list, N, names, title){
    avg_func_n_df(dfs_list, N, names, title, extract_avg_of_n_max_df)
}


get_df_with_max_n_from_grids <- function(n_val, names, complte_dfs_list,average_dfs_list,ward_dfs_list){
    rbind(
        avg_vals_max_n_df(complte_dfs_list,  AVERAGED_N, names,  "complete"),
        avg_vals_max_n_df(average_dfs_list,  AVERAGED_N, names,  "average"),
        avg_vals_max_n_df(ward_dfs_list,     AVERAGED_N, names,  "ward"))
}

# ------------------------------------- Classification Eval - Average of N sampled items from grid

average_N_sample<- function(df, N){
    # Sample N values without repetitions
    sample <- sample(unlist(df), N, replace = FALSE)
    mean(sample, na.rm=TRUE)
}



avg_N_sampled_df <- function(dfs_list, N, names, title){
    avg_func_n_df(dfs_list, N, names, title, average_N_sample)
}


get_df_with_avg_N_sampled_from_grids <- function(n_val, names, complte_dfs_list,average_dfs_list,ward_dfs_list){
    rbind(  avg_N_sampled_df(complte_dfs_list,  AVERAGED_N, names,  "complete"),
        avg_N_sampled_df(average_dfs_list,  AVERAGED_N, names,  "average"),
        avg_N_sampled_df(ward_dfs_list,     AVERAGED_N, names,  "ward"))
}



##################################################### PCI ######################################################

get_PCI <- function(profile_df, method="complete") {
    cleaned_scaled_prf = data.frame()
    cleaned_scaled_prf <- profile_df

    d <- dist(t(cleaned_scaled_prf))
    h_cl <- hclust(d, method)
    cluster_labels <- cutree(h_cl, k = 2)  # Assuming you want to separate into 2 clusters (tumor and normal)

    sample_names <- colnames(cleaned_scaled_prf)
    tumor_samples <- sample_names[grep("T\\.bam$",sample_names)]
    normal_samples <- sample_names[grep("N\\.bam$",sample_names)]
    # normal_indices <- grep("N\\.bam$", names(cluster_labels))
    # tumor_indices <- grep("T\\.bam$", names(cluster_labels))

    concordant_pairs <- sum(cluster_labels[tumor_samples]  == 2) + sum(cluster_labels[normal_samples]  == 1)
    discordant_pairs <-  sum(cluster_labels[normal_samples]  == 2) + sum(cluster_labels[tumor_samples]  == 1)
    total_pairs <- concordant_pairs + discordant_pairs
    PCI <- concordant_pairs / total_pairs
    round(PCI,2)
}

cutoff_and_get_PCI <- function(cutoff, profile_df, method = "complete") {
    cleaned_scaled_prf = data.frame()
    cleaned_scaled_prf <- cutoff_scale_complete_cases(cutoff, profile_df)
    get_PCI(cleaned_scaled_prf, method)
}


cutoff_and_top_sds_get_PCI <- function(sds_count, cutoff, profile_df, method="complete") {
    cleaned_scaled_prf <- data.frame(cutoff_and_get_top_sds(sds_count, cutoff, profile_df))
    if (anyNA(cleaned_scaled_prf[1,1])) {
        return(NA)
    }
    get_PCI(cleaned_scaled_prf,method)
}



get_pci_grid_old <- function(top_sds, cutoffs, df, method="complete") {
        return(get_function_grid_old(top_sds, cutoffs, df, cutoff_and_top_sds_get_PCI, method))
}

get_pci_grid <- function(top_sds, cutoffs, df, method="complete") {
        return(get_function_grid(top_sds, cutoffs, df, get_PCI, method))
}


#----------------------------- PCI Classification eval --------------------------

cutoff_and_top_sds_list_get_PCI <- function(sds_count_list, cutoff, profile_df, method="complete"){
    cutoff_and_top_sds_list_exec_func(sds_count_list, cutoff, profile_df, get_PCI, method)
}

avg_pci_sds_df_new <- function(dfs_list, top_sds, names, cutoff){
    avg_func_res_sds_df_new(dfs_list, top_sds, names, get_PCI, cutoff)
}

# get_pci_grid <- function(top_sds, cutoffs, df, method="complete") {
#         return(get_function_grid_new(top_sds, cutoffs, df, get_PCI, method))
# }

#################################### NMI #################################

library(matrixStats)
library(aricode)
get_NMI <- function(profile_df, method="complete") {
    cleaned_scaled_prf = data.frame()
    cleaned_scaled_prf <- profile_df
    d <- dist(t(cleaned_scaled_prf))
    h_cl <- hclust(d, method)
    cluster_labels <- cutree(h_cl, k = 2)  # Assuming you want to separate into 2 clusters (tumor and normal)
    True_Tumor_Normal_vec <- sapply(names(cluster_labels), function(item) {
        ifelse(grepl("T\\.bam$", item), 'T', 'N')
    })
    nmi <- NMI(cluster_labels,True_Tumor_Normal_vec)
    round(nmi,2)
}

cutoff_and_get_NMI <- function(cutoff, profile_df, method = "complete") {
    cleaned_scaled_prf = data.frame()
    cleaned_scaled_prf <- cutoff_scale_complete_cases(cutoff, profile_df)
    get_NMI(cleaned_scaled_prf, method)
}

cutoff_and_top_sds_get_NMI <- function(sds_count, cutoff, profile_df, method="complete") {
    cleaned_scaled_prf <- data.frame(cutoff_and_get_top_sds(sds_count, cutoff, profile_df))
    if (anyNA(cleaned_scaled_prf[1,1])) {
        return(NA)
    }
    get_NMI(cleaned_scaled_prf,method)
}


get_nmi_grid_old <- function(top_sds, cutoffs, df, method="complete") {
        return(get_function_grid_old(top_sds, cutoffs, df, cutoff_and_top_sds_get_NMI, method))
}

get_nmi_grid <- function(top_sds, cutoffs, df, method="complete") {
        return(get_function_grid(top_sds, cutoffs, df, get_NMI, method))
}


####################################  Silhouette #################################
library(cluster)
get_silhouette_score <- function(df) {
    d <- dist(t(df))
    h_cl <- hclust(d)
    si <- silhouette(cutree(h_cl, k=2), d)
    round(summary(si)$avg.width,2)
}

cutoff_and_get_silhouette_score <- function(cutoff, profile_df){
    cleaned_scaled_prf <- data.frame()
    cleaned_scaled_prf <-  cutoff_scale_complete_cases(cutoff, profile_df)
    return(get_silhouette_score(cleaned_scaled_prf))
}

cutoff_and_top_sds_get_silhouette_score <- function(sds_count, cutoff, profile_df) {
    cleaned_scaled_prf <- data.frame(cutoff_and_get_top_sds(sds_count, cutoff, profile_df))
    if (anyNA(cleaned_scaled_prf[1,1])) {
        return(NA)
    }
    get_silhouette_score(cleaned_scaled_prf)
}


get_silhouette_grid_old <- function(top_sds, cutoffs, df) {
        return(get_function_grid_old(top_sds, cutoffs, df, cutoff_and_top_sds_get_silhouette_score))
}

get_silhouette_grid <- function(top_sds, cutoffs, df) {
        return(get_function_grid(top_sds, cutoffs, df, get_silhouette_score))
}

############################################ T-test ############################################################

add_t_test <- function(df){
    tumor_cols <- colnames(df)[grep("T\\.bam$", colnames(df))]
    normal_cols <- colnames(df)[grep("N\\.bam$", colnames(df))]

    p_values <- apply(df, 1, function(row) {
    ttest_result <- t.test(row[tumor_cols], row[normal_cols])
    ttest_result$p.value
    })
    t_test_df = data.frame(df)
    t_test_df["T-Test"] = p_values
    return(t_test_df)
}


get_avg_ttest_profile <- function(cutoff, profile_df) {
    t_test_df <- data.frame()
    t_test_df <- add_t_test(cutoff_scale_complete_cases(cutoff, profile_df))
    round(mean(t_test_df[,"T-Test"]),2)
}



cutoff_and_top_sds_get_avg_ttest <- function(sds_count, cutoff, profile_df) {
    cleaned_scaled_prf <- data.frame(cutoff_and_get_top_sds(sds_count, cutoff, profile_df))
    if (anyNA(cleaned_scaled_prf[1,1])) {
        return(NA)
    }
    t_test_df <- data.frame()
    t_test_df <- add_t_test(cleaned_scaled_prf)
    round(mean(t_test_df[,"T-Test"]),2)
}

get_avg_ttest_grid_old <- function(top_sds, cutoffs, df) {
        return(get_function_grid_old(top_sds, cutoffs, df, cutoff_and_top_sds_get_avg_ttest))
}


get_avg_ttest_grid <- function(top_sds, cutoffs, df) {
        return(get_function_grid(top_sds, cutoffs, df, get_avg_ttest_profile))
}

################################################################################################33

