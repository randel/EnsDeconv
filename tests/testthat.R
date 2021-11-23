library(testthat)
library(deconv.ensemble)

# testdata not defined
# testdata = deconv.ensemble::testdata
#
# params = get_params("none","none","log","singlecell-rna","test-test",10,"p.value")
#
# res_all =  gen_all_res(count_bulk = bulk_count_ros49,meta_bulk = meta_bulk_ros49,ref_matrix = testdata$ref_list$Nowakowski$count_sc,
#                        meta_ref =  testdata$ref_list$Nowakowski$meta_sc,true_frac = ROS_true,params = params,ncv_input =2,
#                        outpath = "D:/ensemble deconvolution/test/",data_name = "test-test",parallel_comp = T,ncore = 2)
#
# res_all1 =  gen_all_res(count_bulk = bulk_count_ros49,meta_bulk = meta_bulk_ros49,ref_matrix = testdata$ref_list$Nowakowski$count_sc,
#                        meta_ref =  testdata$ref_list$Nowakowski$meta_sc,true_frac = ROS_true,params = params,ncv_input =2,
#                        outpath = "D:/ensemble deconvolution/test/",data_name = "test-test",parallel_comp = F,ncore = 2)
#
#
# test_check("deconv.ensemble")
