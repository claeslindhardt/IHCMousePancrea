#libraries used
library(ggplot2)

#======================================

# T-test based on protein/colour intensity in Stains

#---With standard colour vector

#--------------------------------------


#++++++GZMB exocrine++++++++

#Testing all 2 object glasses with Genzyme B colouring



#Create 3 lists of results for Distal, Proximal part of inslets and one for Control inslets
GZMBexdist <- c(0.1971, 0.1483, 0.2291, 0.1664, 0.3053, 0.2086, 0.1928, 0.1516, 0.2027, 0.1527, 0.0999, 0.1667, 0.1683, 0.2995, 0.3618, 0.2403, 0.1678, 0.2012, 0.177, 0.2685, 0.1529, 0.2052, 0.2628, 0.2625, 0.2347, 0.3289)
GZMBexprox <- c(0.1616, 0.1687, 0.1633, 0.1545, 0.3091, 0.1504, 0.1575, 0.093, 0.2031, 0.1395, 0.0712, 0.194, 0.157, 0.247, 0.3624, 0.169, 0.1633, 0.1599, 0.1764, 0.2171, 0.2029, 0.259, 0.2975, 0.2258, 0.2844, 0.2866)
length(GZMBexdist)
length(GZMBexprox)

#Histograms
hist(GZMBexdist, col='blue')
hist(GZMBexprox, col='green', add=TRUE)

#ttest
ttestexGZMB <- t.test(GZMBexdist, GZMBexprox)

# Print the results of the t-tests
print(ttestexGZMB)

#boxplot
data3 <- data.frame(
  Value = c(GZMBexdist, GZMBexprox),
  Group = factor(rep(c("GZMBexdist", "GZMBexprox"), times = c(length(GZMBexdist), length(GZMBexprox))))
)

ggplot(data3, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot with Individual Data Points", x = "Group", y = "Value")


#======================================

#++++++SMA exocrine++++++++

#Testing all 2 object glasses with SMA colouring

SMAexdist <- c(0.2336, 0.2816, 0.309, 0.301, 0.2811, 0.2484, 0.2929, 0.2772, 0.322, 0.2857, 0.3208, 0.2847, 0.2992, 0.3046, 0.3096, 0.287, 0.2976, 0.3324, 0.3182, 0.3232, 0.3051, 0.3075, 0.2731, 0.3009, 0.3031, 0.2389, 0.2609, 0.2735, 0.2801, 0.2584)
SMAexprox <- c(0.2321, 0.2657, 0.2789, 0.3066, 0.2558, 0.2156, 0.2547, 0.2462, 0.2938, 0.2646, 0.2884, 0.2683, 0.257, 0.2643, 0.2639, 0.2775, 0.2673, 0.2706, 0.267, 0.2813, 0.2921, 0.3063, 0.2777, 0.2676, 0.2928, 0.2697, 0.2325, 0.2546, 0.2816, 0.2264)

length(SMAexdist)
length(SMAexprox)

#Create 3 lists of results for Distal, Proximal part of inslets and one for Control inslets


#Histograms

hist(SMAexprox, col='blue', ylim=c(0,14))

hist(SMAexdist, col='green', add=TRUE)

# Perform 3 t-tests: Dista vs Prox, Dista vs Ctr, Prox vs Ctr
ttestexSMA <- t.test(SMAexprox, SMAexdist)


# Print the results of the t-tests
print(ttestexSMA)

#boxplot
data4 <- data.frame(
  Value = c(SMAexdist, SMAexprox),
  Group = factor(rep(c("SMAexdist", "SMAexprox"), times = c(length(SMAexdist), length(SMAexprox))))
)

ggplot(data4, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot with Individual Data Points", x = "Group", y = "Value")





#++++++SMA islets++++++++

#Testing all 5 object glasses with SMA colouring

#Create 3 lists of results for Distal, Proximal part of inslets and one for Control inslets
SMAisletctr <- c(0.3143, 0.3605, 0.4072, 0.1543, 0.3023, 0.271, 0.4634, 0.2386, 0.3401, 0.3577, 0.213, 0.4408, 0.2058, 0.3949, 0.1901, 0.4532, 0.382, 0.1951, 0.2669, 0.1837, 0.3463, 0.2262, 0.1794, 0.1977, 0.1752, 0.2031, 0.1891, 0.441, 0.1942, 0.2848, 0.3081, 0.1774, 0.4649, 0.1794, 0.477, 0.2037, 0.4738, 0.2125, 0.2506, 0.2472, 0.1661, 0.1945, 0.1766, 0.1846, 0.523)
SMAisletprox <- c(0.4092, 0.3322, 0.2591, 0.3851, 0.547, 0.1416, 0.3662, 0.2803, 0.3325, 0.2965, 0.1608, 0.4136, 0.4385, 0.3089, 0.164, 0.4212, 0.3177, 0.1577, 0.3884, 0.293, 0.2232, 0.3696, 0.2079, 0.4751, 0.1502, 0.2942, 0.2784, 0.2302, 0.1814, 0.2155, 0.1979, 0.1756, 0.1814, 0.4474, 0.1687, 0.2199, 0.1861, 0.4324, 0.1812, 0.5269, 0.1656, 0.5856, 0.1858, 0.1833, 0.1682, 0.2038, 0.2112, 0.2214, 0.1829, 0.1801, 0.1946, 0.1475, 0.2008, 0.2589)
SMAisletdist <- c(0.48, 0.3218, 0.2866, 0.4033, 0.5387, 0.1625, 0.3181, 0.3459, 0.3113, 0.3501, 0.183, 0.3941, 0.4483, 0.3352, 0.1603, 0.4058, 0.1663, 0.4571, 0.3363, 0.243, 0.2969, 0.2293, 0.4929, 0.1679, 0.2969, 0.2798, 0.2206, 0.1979, 0.2112, 0.1885, 0.1738, 0.1735, 0.3723, 0.1949, 0.2401, 0.2142, 0.3621, 0.1874, 0.4135, 0.1776, 0.5086, 0.1879, 0.1858, 0.1847, 0.2142, 0.229, 0.1897, 0.1758, 0.2088, 0.2089, 0.1953, 0.1894, 0.2431, 0.3141)

length(SMAisletctr)
length(SMAisletprox)
length(SMAisletdist)

#Histograms
hist(SMAisletctr, col='blue')
hist(SMAisletprox, col='green', add=TRUE)
hist(SMAisletdist, col='red', add=TRUE)

# Perform 3 t-tests: Dista vs Prox, Dista vs Ctr, Prox vs Ctr
ttest1SMAis <- t.test(SMAisletdist, SMAisletprox)
ttest2SMAis <- t.test(SMAisletdist, SMAisletctr)
ttest3SMAis <- t.test(SMAisletprox, SMAisletctr)

# Print the results of the t-tests
print(ttest1SMAis)
print(ttest2SMAis)
print(ttest3SMAis)

#boxplot
data2 <- data.frame(
  Value = c(SMAisletctr, SMAisletdist, SMAisletprox),
  Group = factor(rep(c("SMAisletctr", "SMAisletdist", "SMAisletprox"), times = c(length(SMAisletctr), length(SMAisletdist), length(SMAisletprox))))
)

ggplot(data2, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot with Individual Data Points", x = "Group", y = "Value")



#======================================


#++++++GZMB islet++++++++

#Testing all 5 object glasses with SMA colouring

#Create 3 lists of results for Distal, Proximal part of inslets and one for Control inslets
GZMBisletctr <- c(0.3781,0.1488,0.2047,0.3303,0.3048,0.2121,0.2604,0.2094,0.3546,0.3329,0.2514,0.2476,0.4457,0.1895,0.2262,0.2754,0.2833,0.391,0.3439,0.1847,0.2822,0.3269,0.2631,0.2534,0.3702,0.2509,0.256,0.3577,0.1128,0.4348,0.2494,0.4195,0.2456,0.2801,0.2894,0.3098,0.291,0.2731,0.2463,0.3019)
GZMBisletdist <- c(0.3282,0.1951,0.4123,0.1527,0.2176,0.2701,0.235,0.3839,0.338,0.192,0.2541,0.2861,0.1706,0.1698,0.2408,0.3134,0.3489,0.2023,0.3349,0.2872,0.2824,0.3082,0.2638,0.2381,0.1215,0.1052,0.2438,0.2832,0.1966,0.3461)
GZMBisletprox <- c(0.2247,0.3359,0.398,0.1555,0.1714,0.2634,0.2282,0.3777,0.3156,0.1811,0.2701,0.2497,0.1798,0.221,0.2614,0.2686,0.3578,0.2029,0.3403,0.2843,0.2733,0.374,0.2673,0.287,0.1705,0.1374,0.2568,0.2871,0.2474,0.4109)


#Histograms
hist(GZMBisletctr, col='blue')

hist(GZMBisletdist, col='green', add=TRUE)

hist(GZMBisletprox, col='red', add=TRUE)

# Perform 3 t-tests: Dista vs Prox, Dista vs Ctr, Prox vs Ctr
ttest1GZMBis <- t.test(GZMBisletdist, GZMBisletprox)
ttest2GZMBis <- t.test(GZMBisletdist, GZMBisletctr)
ttest3GZMBis <- t.test(GZMBisletprox, GZMBisletctr)

# Print the results of the t-tests
print(ttest1GZMBis)
print(ttest2GZMBis)
print(ttest3GZMBis)

#boxplot
# Combine data frames into one
data <- data.frame(
  Value = c(GZMBisletctr, GZMBisletdist, GZMBisletprox),
  Group = factor(rep(c("GZMBisletctr", "GZMBisletdist", "GZMBisletprox"), times = c(length(GZMBisletctr), length(GZMBisletdist), length(GZMBisletprox))))
)

ggplot(data, aes(x = Group, y = Value, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.3) +
  theme_minimal() +
  labs(title = "Boxplot with Individual Data Points", x = "Group", y = "Value")







#Visualisation (boxplots) export
#tiff("ex-CD4.tiff", units="cm", width=6, height=9, res=800, compression = "lzw")
##insert
#dev.off()



