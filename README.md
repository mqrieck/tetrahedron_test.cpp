# tetrahedron_test

Test code for results in a couple papers concerned with the constraints on the angles of a tetrahedron.

Question: Given a fixed triangle in space, when you view its vertices from another point in space, what constraints exist among the angles between the the rays from the viewing point to the vertices?

This simple question has a surprisingly complicated answer. The details for the case when the base triangle is acute have been completely worked out. For an obtuse base triangle, the story is more complicated, but much progress has been made here too.

The programs tetrahedon_test.cpp and dynamic_test_tetrahedron.cpp test my claims made in (yet-to-be-published) research papers dealing with the construction of a tetrahedron based on a given (acute) triangle ABC as a face. The possible interior angles (alpha, beta, gamma) at the new vertex P are restricted in various ways, by a system of (conditional) inequalities. See the "notes" at the top of each source code file, as well as the displayed instructions, for further information.

After compiling either program, one can either execute it without command-line parameters, to indicate that ABC should be an equilateral triangle, or else by specifying three positive integer parameters to indicate the proportions of the interior angles of the triangle ABC. To test an acute base triangle (a completely understood case), be sure that each of these numbers does not exceed the sum of the other two. 

A new source code file, dynamic_test_tetrahedron_obtuse.cpp, has been added. It is tailored to better handle an obtuse base triangle, rather than an acute base triangle. The obtuse base triangle case is significantly more complicated than the acute base triangle case, and this is still a work in progress.

Some animated GIFs are also provided here. Four of these (animated_test_1-4.gif) show examples of executing dynamic_tetrahedon.cpp. One (animated_test_5.gif) shows an example of executing dynamic_tetrahedon_obtuse.cpp. One (spike.gif) shows a 3D image (from Mathematica) of a bounding region (gray) surrounding possible (alpha, beta, gamma) points (purple), for an acute triangle case. Two more (acute.gif and obtuse.gif) show the "Grunert discriminant surface" in "cosines space" (cosines of viewing angles), for quite different base triangle cases. 
