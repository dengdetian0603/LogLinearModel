library(mipfp)

# Example 2: 3-way table (V1,V2,V3) of dim=(2,4,2)
# seed
seed.3d <- array(1,c(2,4,2))
seed.3d[1,1,1] <- 4
seed.3d[1,3,1] <- 10
seed.3d[1,4,2] <- 6
# desired targets (margins) : V1 and (V2,V3)
target.V1 <- c(50, 16)
target.V2.V3 <- array(4, dim=c(4,2))
target.V2.V3[1,1] <- 10
target.V2.V3[3,1] <- 22
target.V2.V3[4,2] <- 14
# list of dimensions of each marginal constrain
tgt.data.3d <- list(target.V1, target.V2.V3)
# storing the description of target data in a list
tgt.list.3d <- list( 1, c(2,3) ,c(1,3))
# calling the Ipfp function
res.3d <- Ipfp(seed.3d, tgt.list.3d, tgt.data.3d, iter=50, print=TRUE, tol=1e-5)
res.3d$xi.hat
loglin(seed.3d,tgt.list.3d,iter=50,param=TRUE)

# ------------------------------------------------------------------------- #
tab0 = ftable(addmargins(HairEyeColor))
mod1 = loglin(HairEyeColor, list(c('Eye', 'Sex'), c('Hair', 'Sex')), 
               fit=TRUE,param = TRUE)
tab1 = ftable(addmargins(mod1$fit))

mod2 <- loglin(HairEyeColor, list(c('Eye', 'Sex'), c('Hair', 'Sex'), 
                                  c('Hair','Eye')), fit=TRUE,param = TRUE)
tab2 = ftable(addmargins(mod2$fit))

mod3 <- loglin(HairEyeColor, list(c('Eye', 'Sex'), c('Hair', 'Sex'), 
                                  c('Hair','Eye'),c('Eye','Sex','Hair')), fit=TRUE,  param = TRUE)
tab3 = ftable(addmargins(mod3$fit))

mod4 = loglin(HairEyeColor, list('Eye', 'Hair', 'Sex'), 
              fit=TRUE,param = TRUE)
tab4 = ftable(addmargins(mod4$fit))



round(tab4-tab0,1)



