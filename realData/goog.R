rm(list=ls(all=TRUE))

# goog.txt has 10 days of data with frequency = 1 minute
# google1 = read.table("goog.txt", sep = ",", header = TRUE)
# head(google1)
# dim(google1)
# len = length(google1[[1]])
# timemat = numeric(length = len - 1)

# for(i in c(1:len))
# {
# 	timemat[i] = google1[i+1,1] - google1[i,1] 
# }

# hist(timemat)

# goog2.txt has 25 days of data with frequency = 1 minute
# google2 = read.table("goog2.txt", sep = ",", header = TRUE)
# google2 = google2[1:100,]
# head(google2)
# dim(google2)
# len = length(google2[[1]])
# timemat = numeric(length = len - 1)

# for(i in c(1:len))
# {
# 	timemat[i] = google2[i+1,1] - google2[i,1] 
# }

# plot(timemat)

# goog3.txt has 25 days of data with frequency = 1 hour
# google3 = read.table("goog3.txt", sep = ",", header = TRUE)
# # google3 = google3[1:20,]
# head(google3)
# dim(google3)
# len = length(google3[[1]])
# timemat = numeric(length = len - 1)

# for(i in c(1:len))
# {
# 	timemat[i] = google3[i+1,1] - google3[i,1] 
# }

# plot(timemat)

# goog4.txt has 100 days of data with frequency = 600 seconds
# Total minutes in a trading day = 390

google4 = read.table("goog4.txt", sep = ",", header = TRUE)
head(google4)
dim(google4)
len = length(google4[[1]])
timemat = numeric(length = len - 1)

for(i in c(1:len))
{
	timemat[i] = google4[i+1,1] - google4[i,1] 
}

# plot(timemat)
google4 = google4[c(-3, -4, -6)]

timestamp = 40
days = dim(google4)[1] / timestamp
# one matrix which has closing stock prices and the time for 100 days
closingprice = matrix(0, nrow = days, ncol = timestamp)
starttime = 1
stoptime = timestamp

for(i in c(1:days))
{
	closingprice[i,] = google4[starttime:stoptime,2]
	starttime = stoptime + 1
	stoptime = (i+1) * timestamp
}

openingprice = matrix(0, nrow = days, ncol = timestamp)
starttime = 1
stoptime = timestamp

for(i in c(1:days))
{
	openingprice[i,] = google4[starttime:stoptime,3]
	starttime = stoptime + 1
	stoptime = (i+1) * timestamp
}

# mylist = list(openingprice, closingprice)
save(closingprice, file = 'goog4.RData')