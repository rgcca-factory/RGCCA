
#' European Customer Satisfaction Index
#'
#' @format A data frame with 250 rows and 24 variables
#' @docType data
#'
#' @description
#' The European Consumer Satisfaction Index (ECSI) is an economic indicator that
#' measures customer satisfaction. ECSI is an adaptation of the Swedish Customer
#' Satisfaction Barometer (Fornell, 1992) and is compatible with the American
#' Customer Satisfaction Index.  The indicators describing the latent variables
#' are given for the Mobile Phone Industry. The original items scaled from 1 to 10
#' have been transformed into new normalized variables. The minimum possible value
#' of each variable is 0 and its maximum possible value is equal to 10.
#'
#' \describe{
#' - IMAG: Image of the phone provider (eta_1) \cr
#' \verb{    }(a) reputation of the phone provider \cr
#' \verb{    }(b) trustworthiness \cr
#' \verb{    }(c) seriousness \cr
#' \verb{    }(d) solidness,  \cr
#' \verb{    }(e) caring about customer's needs.\cr
#'
#' - EXPE: Customer Expectations of the overall quality (eta_2)\cr
#' \verb{    }(a)   Expectations for the overall quality of your "mobile phone
#' provider" at the moment you became customer of this provider \cr
#' \verb{    }(b)   Expectations for your "mobile phone provider" to provide
#' products and services to meet your personal need \cr
#' \verb{    }(c)   How often did you expect that things could go wrong at your
#' "mobile phone provider" \cr
#'
#' - QUAL: Perceived Quality (eta_3) \cr
#' \verb{    }(a) Overall perceived quality \cr
#' \verb{    }(b) Technical quality of the network \cr
#' \verb{    }(c) Customer service and personal advice offered \cr
#' \verb{    }(d) Quality of the services you use \cr
#' \verb{    }(e) Range of services and products offered \cr
#' \verb{    }(f) Reliability and accuracy of the products and services provided \cr
#' \verb{    }(g) Clarity and transparency of information provided \cr
#'
#' - VAL: Perceived Value (eta_4)\cr
#' \verb{    }(a)   Given the quality of the products and services offered by
#' your "mobile phone provider" how would you rate the fees and prices that
#' you pay for them?\cr
#' \verb{    }(b)   Given the fees and prices that you pay for your mobile phone
#' provider how would you rate the quality of the products and services offered
#' by your "mobile phone provider"? \cr
#'
#' - SAT: Customer Satisfaction (eta_5)\cr
#' \verb{    }(a)   Overall satisfaction \cr
#' \verb{    }(b)   Fulfillment of expectations \cr
#' \verb{    }(c)   How well do you think your "mobile phone provider" compares
#' with your ideal "mobile phone provider"? \cr
#'
#' - LOY: Customer Loyalty (eta_6)   \cr
#' \verb{    }(a)   If you would need to choose a new "mobile phone provider" how
#' likely is it that you would choose your provider again? \cr
#' \verb{    }(b)   Let us now suppose that other "mobile phone provider"s decide to
#' lower their fees and prices, but your "mobile phone provider" stays at the same
#' level as today. At which level of difference (in \%) would you choose another
#' "mobile phone provider"? \cr
#' \verb{    }(c)   If a friend or colleague asks you for advice, how likely is it
#' that you would recommend your "mobile phone provider"? }
#'
#' @references Fornell C. (1992): A national customer satisfaction barometer.
#' The Swedish experience. Journal of Marketing, (56), 6-21.
#' @usage data(ECSI)
#' @keywords datasets
"ECSI"
