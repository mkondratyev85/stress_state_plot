# stress_state_plot

Example of usage:

```
python main.py --sigma1_orientation_and_sigma3_direction="125 34 322" --pressure=-40 --mu_sigma=0 --tau=1 --png_path=/tmp/fig.png
```

It is possible to plot set of fractures against given stress state:

```
python main.py --sigma1_orientation_and_sigma3_direction="125 34 322" --pressure=-40 --mu_sigma=0 --tau=1 --fractures=fractures.txt --png_path=/tmp/fig.png
```

<img src="figs/fig.png" width="640" />

<img src="figs/fig_morh.jpg" width="640" />

Example of usage in interactive mode:

```
python main.py --sigma1_orientation_and_sigma3_direction="125 34 322" --pressure=-40 --mu_sigma=0 --tau=1 --fractures=fractures.txt --gui
```

<img src="figs/interactive_screencast.webm" width="640" />

