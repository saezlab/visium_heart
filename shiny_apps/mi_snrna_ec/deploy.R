library('rsconnect')

tkn <- readRDS(file = "./tkn.rds")

rsconnect::setAccountInfo(name = tkn$name,
                          token= tkn$token,
                          secret= tkn$secret)

options(repos = BiocManager::repositories())
options(shiny.maxRequestSize=500*1024^2)

rsconnect::deployApp(appFiles = c('app.R', "testcells.rds"),
                     appName = 'test_app', account = 'saezlab', server = 'shinyapps.io')
