# OnePassStatistics
A one pass statistics class for MATLAB

General usage is very simple, just once initialize:

```stat = OnePassStatistics;```

Within the loop, add the current data iterate:

```stat = stat.addData(data);```

And whenever you want ask for the current sample statistics:
```
stat.mean();
stat.var();
stat.std();
stat.cov();
stat.skewness();
stat.kurtosis();
```

Have Fun!

To Do: Weighted statistics!
