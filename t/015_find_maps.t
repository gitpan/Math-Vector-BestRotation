#!perl -T

use strict;
use warnings;

use Test::More tests => 51;
use Math::MatrixReal;

BEGIN {
    use_ok('Math::Vector::BestRotation');
}

sub get_vector_pair {
    my ($matrix, $coords) = @_;
    my $mapped;

    $coords = Math::MatrixReal->new_from_cols([$coords])
	if(ref($coords) eq 'ARRAY');
    $mapped = $matrix * $coords;
    return([map { $coords->element($_, 1) } (1, 2, 3)],
	   [map { $mapped->element($_, 1) } (1, 2, 3)]);
}

sub best_orthogonal {
    my $rot = Math::Vector::BestRotation->new;
    my $matrix;
    my $ref;
    my $diff;
    my ($rot1, $rot2, $rot3);

    can_ok($rot, 'best_orthogonal');

    my $flip = Math::MatrixReal->new_from_string(<<'    MATRIX');
        [ -1   0   0  ]
	[  0  -1   0  ]
	[  0   0  -1  ]
    MATRIX

    $ref = Math::MatrixReal->new_from_string(<<'    MATRIX');
        [  0  -1   0  ]
	[  1   0   0  ]
	[  0   0   1  ]
    MATRIX
    $rot->add_pair([1, 0, 0], [0, 1, 0]);
    $rot->add_pair([0, 1, 0], [-1, 0, 0]);
    $matrix = $rot->best_orthogonal;
    isa_ok($matrix, 'Math::MatrixReal');
    cmp_ok(abs($matrix->det - 1), '<', 1e-9, 'is special orthogonal');
    $diff = $matrix - $ref;
    for(my $i=0;$i<3;$i++) {
	for(my $j=0;$j<3;$j++) {
	    cmp_ok(abs($diff->element($i+1, $j+1)), '<', 1e-9,
		   sprintf("result from scratch %d %d", $i+1, $j+1));
	}
    }

    # another simple one
    $rot->clear;
    $ref = Math::MatrixReal->new_from_string(<<'    MATRIX');
        [  1   0   0  ]
	[  0   0   1  ]
	[  0  -1   0  ]
    MATRIX
    $rot->add_pair(get_vector_pair($ref, [0, 1, 0]));
    $rot->add_pair(get_vector_pair($ref, [0, 0, 1]));
    $matrix = $rot->best_orthogonal;
    isa_ok($matrix, 'Math::MatrixReal');
    cmp_ok(abs($matrix->det - 1), '<', 1e-9, 'is special orthogonal');
    $diff = $matrix - $ref;
    for(my $i=0;$i<3;$i++) {
	for(my $j=0;$j<3;$j++) {
	    cmp_ok(abs($diff->element($i+1, $j+1)), '<', 1e-9,
		   sprintf("another simple %d %d", $i+1, $j+1));
	}
    }

    # arbitrary rotation
    $rot->clear;
    $rot1 = Math::MatrixReal->new_from_rows([
        [  1,          0,           0  ],
	[  0,   cos(0.3),   -sin(0.3)  ],
	[  0,   sin(0.3),    cos(0.3)  ],
					    ]);
    $rot2 = Math::MatrixReal->new_from_rows([
	[ cos(2.1),   0,    -sin(2.1)  ],
        [        0,   1,            0  ],
	[ sin(2.1),   0,     cos(2.1)  ],
					    ]);
    $rot3 = Math::MatrixReal->new_from_rows([
	[  cos(-1),   -sin(-1),  0   ],
	[  sin(-1),    cos(-1),  0   ],
        [        0,          0,  1   ],
					    ]);
    $ref = $rot1 * $rot2 * $rot3;
    cmp_ok(abs($ref->det - 1), '<', 1e-9, 'ref is special orthogonal');
    $rot->add_pair(get_vector_pair($ref, [1, 2, 3]));
    $rot->add_pair(get_vector_pair($ref, [4, -1, 0]));
    $rot->add_pair(get_vector_pair($ref, [3, 5, 1]));
    $rot->add_pair(get_vector_pair($ref, [-7, -0.123, -23.786]));
    $matrix = $rot->best_orthogonal;
    isa_ok($matrix, 'Math::MatrixReal');
    cmp_ok(abs($matrix->det - 1), '<', 1e-9, 'is special orthogonal');
    $diff = $matrix - $ref;
    for(my $i=0;$i<3;$i++) {
	for(my $j=0;$j<3;$j++) {
	    cmp_ok(abs($diff->element($i+1, $j+1)), '<', 1e-9,
		   sprintf("arbitrary with 4 pairs %d %d", $i+1, $j+1));
	}
    }

    # and flipped
    $rot->clear;
    $ref = $ref * $flip;
    cmp_ok(abs($ref->det + 1), '<', 1e-9, 'ref has det -1');
    $rot->add_pair(get_vector_pair($ref, [1, 2, 3]));
    $rot->add_pair(get_vector_pair($ref, [4, -1, 0]));
    $rot->add_pair(get_vector_pair($ref, [3, 5, 1]));
    $rot->add_pair(get_vector_pair($ref, [-7, -0.123, -23.786]));
    cmp_ok($rot->matrix_r->det, '<', 0, 'R has det < 0');
    $matrix = $rot->best_orthogonal;
    isa_ok($matrix, 'Math::MatrixReal');
    cmp_ok(abs($matrix->det + 1), '<', 1e-9, 'has det -1');
    $diff = $matrix - $ref;
    for(my $i=0;$i<3;$i++) {
	for(my $j=0;$j<3;$j++) {
	    cmp_ok(abs($diff->element($i+1, $j+1)), '<', 1e-9,
		   sprintf("same one flipped %d %d", $i+1, $j+1));
	}
    }
    $matrix = $rot->best_rotation;
    isa_ok($matrix, 'Math::MatrixReal');
    cmp_ok(abs($matrix->det - 1), '<', 1e-9,
	   'best rotation is special orthogonal');
}

best_orthogonal;

