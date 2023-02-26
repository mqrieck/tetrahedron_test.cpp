# tetrahedron_test

Test code for results in a couple papers concerned with the constraints on the angles of a tetrahedron. 

THIS MATERIAL IS INTIMATELY RELATED TO THE P3P PROBLEM.

**Question:** Given a fixed triangle in space, when you view its vertices from another point in space, what constraints exist among the angles (α, β and γ) between the rays from the viewing point to the vertices?

This simple question has a surprisingly complicated, yet still tractable (in a sense), answer. The details for the case when the base triangle is acute have been completely worked out. Acute_Triangle_Case.pdf is a fact sheet for this case. For an obtuse base triangle, the story is more complicated, but much progress has been made here too.

**Recent news:** There is now a link to a YouTube video that discusses this and a related problem, and a link to the Mathematica notebook used in that video, and a "four problems" fact sheet covering the same material. OL1.gif and OL2.gif show off some subtle aspects of the OL Problem, without really explaining these though. ObtuseK.gif is quite a nice 3D animation showing the surface of which the discriminant vanishes in cosines space (also not explained here). 

**Older material:**

The programs tetrahedon_test.cpp and dynamic_test_tetrahedron.cpp test the claims in the fact sheet and in (yet-to-be-published) research papers dealing with the construction of a tetrahedron based on a given acute triangle ABC as a face. The possible interior angles (α, β and γ) at the new vertex P are restricted in various ways, by a system of (conditional) inequalities. For more information, see the fact sheet, the notes at the top of each source code file, and the displayed instructions when the programs are run.

After compiling either program, one can either execute it without command-line parameters, to indicate that ABC should be an equilateral triangle, or else specify three positive integer command-line parameters to indicate the proportions of the interior angles of the triangle ABC. To test an acute base triangle, be sure that each of these numbers does not exceed the sum of the other two. 

More recent source code file, dynamic_test_tetrahedron_obtuse.cpp, has been added. This program is tailored to better handle an obtuse base triangle, rather than an acute base triangle. The obtuse base triangle case is significantly more complicated than the acute base triangle case, and this is still a work in progress. When executing this program, the three command-line numbers are required, and the first of these must be greater than the sum of the other two. 

An animated GIF (rays_problem.gif) shows even more recent experimental evidence for the bounding problem. It show results that are very similar to those of the "dynamic" C++ progam, but Mathematica was used here instead. However, by generating data points in a completely different way, much clearer results were produced. Here we see possible pairs (cos α, cos β) (the blue region) for a given value of γ. Another animation (lines_problem.gif) shows similar experimental results for a related but simpler problem (see below). 

Some other animated GIFs are also provided. Four of these GIFs (animated_test__1-4.gif) show examples of executing dynamic_tetrahedon.cpp. Two GIFs (animated_test__5-6.gif) shows examples of executing dynamic_tetrahedon_obtuse.cpp. One GIF (spike.gif) shows a 3D image of the bounding region (gray) given by the inequalities in the fact sheet, surrounding some possible (α,β,γ) points (purple), in "view angle space" (coordinates being the "view angles" α, β and γ), for an acute triangle case. Two more GIFs (acute.gif and obtuse.gif) show the "discriminant surface" of "Grunert's system of equations," in "cosines space" (the coordinates being the three cosines of the view angles), the first animation being for an acute triangle, and the second being for an obtuse triangle. 

Another GIF (overlay.gif) shows that by taking four copies of a "spike" (as in spike.gif), these being rotations of each other, they combine to fill out a regular tetrahedron surface, but generally not its interior. The resulting solid region corresponds to a related, but simpler problem (the "lines problem"), i.e. the problem of trying to fit a given triangle so that its vertices touch corresponding fixed lines (rather than rays) through a fixed point. 
