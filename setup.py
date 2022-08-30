from setuptools import find_namespace_packages, setup

setup(name='stress_state_plot',
      version='0.1',
      packages=find_namespace_packages(where="src"),
      package_dir={"": "src"},
      package_data={
          "canvas.patterns": ["*.txt"],
          },
      install_requires=["click",
                        "scipy",
                        "numpy",
                        "matplotlib",
                        ],
      entry_points={
        'console_scripts': ['stress_state_plot=stress_state_plot.__main__:run', ]
      },
      )
