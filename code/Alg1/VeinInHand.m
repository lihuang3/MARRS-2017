figure(2);clf
axis([0 220 0 140]);

x_obs = [ 0 16 48 69  0 0  NaN 71 95 80 71 NaN  0 30  0  0 0 NaN  0 40 90 60 0 0 NaN...
     0   60 110 140 0  0   NaN 30 42 100 106.5 75 30 NaN 61 69 100 95 128 116 61 NaN ...
    87 100 150 159 116 87 NaN 132 157 140 132 NaN ...
    88 113.5 147  128 118 88 NaN 145 200 200 170 162 185 210 214 196 210 201 186 119 ...
     160  157 220   220 157 145 NaN 157 220 220 198 195 162 127 157 NaN ...
     73 79 162 181 181 140 145  130 146 145  115 73 NaN 173 189 181 181 220 220 198 173];

y_obs = [50 56 43  0  0 50 NaN 17   0  0 17 NaN 60 70 80 65 60 NaN 90 75 75 100 100 90 NaN ...
    110 110 124 140 140  110 NaN 61 66 66    60   43 61 NaN 38 29  8  0   0  15 38 NaN...
    38 33.5 39  54  54 38 NaN  10   0   0  10 NaN...
    90  94.5  75   64  64 90 NaN  64 64  54  54  40  43  54  45  37   20  15  33 26 ...
    10  0   0    75  70  64     NaN 80  85 104  104  95 103  97 80 NaN ...
    103 98 113 138 140 140 132  124 124 114 114 103 NaN 111 132 138 140 140 104 104 111];
    
mapshow(x_obs,y_obs,'DisplayType','polygon',...
    'FaceColor',[182,228,255]./255,'LineStyle','none')
  