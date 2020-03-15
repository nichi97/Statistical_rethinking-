output_list <- list()
for(i in seq(1,10)){
    d <- runif(40,0,2)
    output_list[[i]] <- d
}

mean_total <- list()
for(l in seq(1, length(output_list))){
  mean_list <- vector("numeric")
  curr_ls <- output_list[[l]]
  for(i in seq(length(curr_ls))){
    mean_list[[i]] <- mean(curr_ls[1:i])
  }
  mean_total[[l]] <- mean_list
}


  