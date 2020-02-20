 s = '1,2, 3'
res= char_to_list(s)
 test_that("test_char_to_list",{
     expect_true(res[2]=="2")
 })
 