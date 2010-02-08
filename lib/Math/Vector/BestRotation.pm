package Math::Vector::BestRotation;

use warnings;
use strict;

use 5.008008;

use Carp;
use List::Util qw(sum);
use Math::MatrixReal;

=head1 NAME

C<Math::Vector::BestRotation> - best rotation to match two vector sets

=head1 VERSION

Version 0.007

=cut

our $VERSION = '0.007';


###########################################################################
#                                                                         #
#                             Init Process                                #
#                                                                         #
###########################################################################

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;
    $self->init(@args);
    return $self;
}

sub init {
    my ($self, %args) = @_;

    $self->{matrix_r} = [0, 0, 0, 0, 0, 0, 0, 0, 0];

    foreach(keys %args) {
	my $meth = $_;
	if($self->can($meth)) { $self->$meth($args{$meth}) }
	else { carp "Unrecognized init parameter $meth.\n" }
    }
}

###########################################################################
#                                                                         #
#                               Accessors                                 # 
#                                                                         #
###########################################################################

sub matrix_r {
    my ($self, @args) = @_;

    croak "Attribute matrix_r is readonly.\n" if(@args);
    return Math::MatrixReal->new_from_rows
	([[$self->{matrix_r}->[0],
	   $self->{matrix_r}->[1],
	   $self->{matrix_r}->[2]], 
	  [$self->{matrix_r}->[3],
	   $self->{matrix_r}->[4],
	   $self->{matrix_r}->[5]],
	  [$self->{matrix_r}->[6],
	   $self->{matrix_r}->[7],
	   $self->{matrix_r}->[8]]]);
}

###########################################################################
#                                                                         #
#                             Retrieve Data                               #
#                                                                         #
###########################################################################

sub _compute_bases {
    my ($self)   = @_;
    my $r        = $self->matrix_r;
    my $rtr      = (~$r) * $r;
    my ($mu, $a) = $rtr->sym_diagonalize;

    my $sorted = 0;
    while(!$sorted) {
	$sorted = 1;
	if(abs($mu->element(1, 1)) < abs($mu->element(2, 1))) {
	    $sorted = 0;
	    $mu->swap_row(1, 2);
	    $a->swap_col(1, 2);
	}
	if(abs($mu->element(2, 1)) < abs($mu->element(3, 1))) {
	    $sorted = 0;
	    $mu->swap_row(2, 3);
	    $a->swap_col(2, 3);
	}
    }

    # make sure, eigensystem is right-handed
    my $a3 = $a->column(1)->vector_product($a->column(2));
    $a->assign(1, 3, $a3->element(1, 1));
    $a->assign(2, 3, $a3->element(2, 1));
    $a->assign(3, 3, $a3->element(3, 1));

    if($mu->element(1, 1) < 1e-12) {
	carp("The largest eigenvalue of R^t*R is smaller than ".
	     "1e-12. The code will try to go on, but may die with ".
	     "division by zero later. In any case, numerical ".
	     "instability of the result has to be expected.\n");
    }

    my $b = [];  # matrix of column vectors b_k
    my $v;       # vector buffer
    my $l;       # length buffer
    my $col;     # column matrix buffer

    # determine b_1
    $col = $r * $a->column(1);
    $v   = [map { $col->element($_, 1) } (1, 2, 3)];
    $l   = sqrt(sum(map { $v->[$_]**2 } (0, 1, 2)));

    # If length is zero all b's are arbitrary. Otherwise normalize.
    if($l == 0) {
	carp("R is (too close to) the nullmatrix. ".
	     "The result will be arbitrary.\n");
	$v = [1, 0, 0];
    }
    else { @$v = map { $_ / $l } @$v }
    @$b[0, 3, 6] = @$v;

    # determine b_2
    $col = $r * $a->column(2);
    $v   = [map { $col->element($_, 1) } (1, 2, 3)];
    $l   = sqrt(sum(map { $v->[$_]**2 } (0, 1, 2)));

    # if length is zero, we choose a vector perpendicular to b_1
    if($l == 0) {
	carp("The two smaller eigenvalues are 0. The result is not ".
	     "unique.\n");

	if(abs($b->[0]) < abs($b->[3])) {
	    if(abs($b->[0]) < abs($b->[6])) {
		$b->[1] = 0;
		$b->[4] = -$b->[6];
		$b->[7] = $b->[3];
	    }
	    else {
		$b->[1] = -$b->[3];
		$b->[4] = $b->[0];
		$b->[7] = 0;
	    }
	}
	else {
	    if(abs($b->[3]) < abs($b->[6])) {
		$b->[1] = -$b->[6];
		$b->[4] = 0;
		$b->[7] = $b->[0];
	    }
	    else {
		$b->[1] = -$b->[3];
		$b->[4] = $b->[0];
		$b->[7] = 0;
	    }
	}
    }
    # otherwise we normalize carefully
    else {
	@$v = map { $_ / $l } @$v;

	# $l could have been very small, though.
	# Therefore, we make b_2 orthogonal to b_1 and normalize again.
	my $p = $v->[0] * $b->[0] + $v->[1] * $b->[3] + $v->[2] * $b->[6];
	@$v = map { $v->[$_] - $p * $b->[3*$_] } (0, 1, 2);
	$l  = sqrt(sum(map { $v->[$_]**2 } (0, 1, 2)));
	@$b[1, 4, 7] = map { $_ / $l } @$v;
    }

    # determine b_3 as by cross product
    $b->[2] = $b->[3] * $b->[7] - $b->[6] * $b->[4];
    $b->[5] = $b->[6] * $b->[1] - $b->[0] * $b->[7];
    $b->[8] = $b->[0] * $b->[4] - $b->[3] * $b->[1];

    return($a, Math::MatrixReal->new_from_rows
	([[@$b[0, 1, 2]], [@$b[3, 4, 5]], [@$b[6, 7, 8]]]));
}

sub _matrix_a_to_b {
    my ($self, $a, $b) = @_;

    return($b * (~$a));
}

sub best_orthogonal {
    my ($self)  = @_;
    my ($a, $b) = $self->_compute_bases;

    if($self->matrix_r->det < 0) {
	$b = $b * Math::MatrixReal->new_from_rows
	    ([[1, 0, 0], [0, 1, 0], [0, 0, -1]]);
    }
    return $self->_matrix_a_to_b($a, $b);
}

sub best_rotation {
    my ($self)  = @_;
    my ($a, $b) = $self->_compute_bases;

    return $self->_matrix_a_to_b($a, $b);
}

sub best_flipped_rotation {
    my ($self)  = @_;
    my ($a, $b) = $self->_compute_bases;

    $b = $b * Math::MatrixReal->new_from_rows
	([[1, 0, 0], [0, 1, 0], [0, 0, -1]]);
    return $self->_matrix_a_to_b($a, $b);
}

###########################################################################
#                                                                         #
#                              Change Data                                # 
#                                                                         #
###########################################################################

sub add_pair {
    my ($self, $x, $y) = @_;
    my $matrix_r       = $self->{matrix_r};

    for(my $i=0;$i<3;$i++) {
	for(my $j=0;$j<3;$j++) {
	    $matrix_r->[3*$i+$j] += $y->[$i] * $x->[$j];
	}
    }
}

sub add_many_pairs {
    my ($self, $x, $y) = @_;
    my $matrix_r       = $self->{matrix_r};

    for(my $n=0;$n<@$x;$n++) {
	for(my $i=0;$i<3;$i++) {
	    for(my $j=0;$j<3;$j++) {
		$matrix_r->[3*$i+$j] += $y->[$n]->[$i] * $x->[$n]->[$j];
	    }
	}
    }
}

sub clear {
    my ($self) = @_;

    $self->{matrix_r} = [0, 0, 0, 0, 0, 0, 0, 0, 0];
}


1;

__END__

=head1 SYNOPSIS

    use Math::Vector::BestRotation;

    my $best = Math::Vector::BestRotation->new();

    $best->add_pair([1, 2, 3], [4, 5, 6]);
    $best->add_pair([0, -1, 5], [-3, 6, 0]);
    .
    .
    .
    $best->add_pair([3, -1, 5], [0.1, -0.7, 4]);

    my $ortho = $best->best_orthogonal;
    my $rot   = $best->best_rotation;
    my $flip  = $best->best_flipped_rotation;

    # start over
    $best->clear;

=head1 DESCRIPTION

Assume that you have a list of vectors v_1, v_2, v_3, ..., v_n and
an equally sized list of vectors w_1, w_2, ..., w_n. A way to
quantify how similar these lists are to each other is to compute the
sum of the squared distances between the vectors:
sum((w_1 - v_1)**2 + ... + (w_n - v_n)**2).
In the literature, this sum is sometimes divided by 2 or divided
by n or divided by n and the square root is taken ("root mean
square" or RMS deviation).

In some situations, one data set can be arbitrarily rotated with
respect to the other one. In this case, one of them has to be
rotated in order to calculate the RMS deviation in a meaningful way.
C<Math::Vector::BestRotation> solves this problem. It calculates the
best orthogonal map U between the v_i and w_i. "Best" means here
that the RMS deviation between Uv and w as calculated above is
minimized.

An orthogonal map can be a rotation or a rotation combined with
a reflection (I call that a flipped rotation). This module enables
you to find the best orthogonal map, the best rotation, or the
best flipped rotation between two given vector sets.

The algorithm implemented here is based on two research papers
listed in the L<ACKNOLEDGEMENT|/ACKNOWLEDGEMENT> section. It works
for higher dimensional vector spaces as well, but the current
implementation supports only three-dimensional vectors. This
limitation is going to be remedied in a future version of this
module.

The two data sets could not only rotated with respect to each
other, but also translated. This translation can be removed prior
to the determination of the rotation by aligning the centers of
masses of the two vector sets. However, this procedure is not
offered by C<Math::Vector::BestRotation> and possibly will never
be, because this would require to store the full data sets in
memory which is not necessary now.

=head1 INTERFACE

=head2 Constructors

=head3 new

  Usage   : Math::Vector::BestRotation->new(%args)
  Function: creates a new Math::Vector::BestRotation object
  Returns : a Math::Vector::BestRotation object
  Args    : initial attribute values as named parameters

Creates a new C<Math::Vector::BestRotation> object and calls
C<init(%args)>. If you subclass C<Math::Vector::BestRotation>
overload C<init>, not C<new>.

=head3 init

  Usage   : only called by new
  Function: initializes attributes
  Returns : nothing
  Args    : initial attribute values as named parameters

If you overload C<init>, your method should also call this one.
It provides the following functions:

=over 4

=item * For each given argument (which has not been deleted by
the previous actions) it calls the accessor with the same name to
initialize the attribute. If such an accessor does not exist a
warning is printed and the argument is ignored.

=back

=head2 Public Attributes

=head3 matrix_r

This attributes holds a matrix built from the pairs of vectors and
used to compute the desired orthogonal map. It is called R in the
documentation and the underlying L<Research papers|/ACKNOWLEDGEMENTS>.
The accessor is readonly. At startup, C<matrix_r> is initialized with
zeros.

Note that the matrix is stored internally as an array to speed up
data acquisition. When you call the accessor a C<Math::MatrixReal>
object is created. This implies that such an object is not updated
as you add more vector pairs. You have to call the accessor again
to get a new object. Accordingly, changing of your retrieved matrix
does not alter the underlying matrix stored in the
C<Math::Vector::BestRotation> object.

=head2 Methods for Users

=head3 add_pair

  Usage   : $obj->add_pair([1, 2, 3], [0, 7, 5])
  Function: updates matrix_r
  Returns : nothing
  Args    : a pair of vectors, each as an ARRAY reference

The orthogonal map computed by this module will try to map the first
vector of each pair onto the corresponding second vector. This method
uses the new vector pair to
update the matrix R which is later used to compute the best map.
The vectors are discarded afterwards and can therefore not be
removed once they have been added.

In some applications, very many vector pairs will be added making
this the rate limiting step of the calculation. Therefore, no
convenience functionality is offered. For example, the method
strictly requires ARRAY references. If your vectors are stored e.g.
as L<Math::VectorReal|Math::VectorReal> objects you have to turn
them into ARRAY references yourself. Furthermore, no error checking
whatsoever is performed. The method expects you to provide valid
data. See also L<add_many_pairs|/add_many_pairs>.

=head3 add_many_pairs

  Usage   : $obj->add_many_pairs([[1, 2, 3], [3, 0, 6]],
                                 [[0, 7, 5], [-2, 1, -1]])
  Function: updates matrix_r
  Returns : nothing
  Args    : a pair of vectors, each as an ARRAY reference

An alternative to L<add_pair|/add_pair>. It expects two lists of
vectors. The first one contains the first vector of each pair,
the second one contains the second vector of each pair (see
L<add_pair|/add_pair>. If you have many vector pairs to add it
is probably faster to build these lists and then use this method
since it saves you a lot of method calls.

For perfomance reasons, no checks are performed not even if the
two lists have equal sizes. You are expected to provide valid
data.

=head3 best_orthogonal

  Usage   : $matrix = $obj->best_orthogonal
  Function: computes the best orthogonal map between the vector sets
  Returns : a Math::MatrixReal object
  Args    : none

Computes the best orthogonal map between the two vector sets, i.e.
the orthogonal map that minimizes the sum of the squared distances
between the image of the first vector of each pair and the
corresponding second vector. This map can be either a rotation or
a rotation followed by a reflection.

The representing matrix of the found map is returned in form of a
L<Math::MatrixReal|Math::MatrixReal> object.

=head3 best_rotation

  Usage   : $matrix = $obj->best_rotation
  Function: computes the best rotation between the vector sets
  Returns : a Math::MatrixReal object
  Args    : none

This is identical to L<best_orthogonal|/best_orthogonal> except that
it finds the best special orthogonal map (this means that the
determinant is +1, i.e. the map is a true rotation).

The method computes the best rotation between the two vector sets,
i.e. the rotation that minimizes the sum of the squared distances
between the image of the first vector of each pair and the
corresponding second vector.

The representing matrix of the found map is returned in form of a
L<Math::MatrixReal|Math::MatrixReal> object.

=head3 best_flipped_rotation

  Usage   : $matrix = $obj->best_flipped_rotation
  Function: computes the best rotation combined with a reflection
  Returns : a Math::MatrixReal object
  Args    : none

This is identical to L<best_orthogonal|/best_orthogonal> except that
it finds the best orthogonal map with determinant is -1. I do not
know why one would want that, but the method is included for
completeness.

The representing matrix of the found map is returned in form of a
L<Math::MatrixReal|Math::MatrixReal> object.

=head3 clear

  Usage   : $obj->clear
  Function: resets the object
  Returns : nothing
  Args    : none

This method resets L<matrix_r|/matrix_r> to the null matrix (all
entries equal zero). This enables you to start from scratch with
two new vector sets without destroying the object.

=head1 DIAGNOSTICS

=head2 Exceptions

Sorry, not documented, yet.

=head2 Warnings

Sorry, not documented, yet.


=head1 BUGS AND LIMITATIONS

=head2 Bugs

No bugs have been reported, but the code should be considered as
beta quality.

Please report any bugs or feature requests to
C<bug-math-vector-bestrotation at rt.cpan.org>, or through
the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math-Vector-BestRotation>.
I will be notified, and then you
will automatically be notified of progress on your bug as I make
changes.

=head2 Limitations

See L<DESCRIPTION|/DESCRIPTION>.


=head1 AUTHOR

Lutz Gehlen, C<< <perl at lutzgehlen.de> >>


=head1 ACKNOWLEDGEMENTS

The algorithm implemented here is based on two research papers by
Wolfgang Kabsch, Max-Planck-Institut für Medizinische Forschung,
Heidelberg, Germany:

=over 4

=item [1] Kabsch, W. (1976). A solution for the best rotation to
relate two sets of vectors. Acta Cryst., A32, 922

=item [2] Kabsch, W. (1978). A discussion of the solution for the
best rotation to relate two sets of vectors. Acta Cryst., A34,
827-828

=back

=head1 LICENSE AND COPYRIGHT

Copyright 2010 Lutz Gehlen.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=cut
