
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<p>
Several R packages providing (algorithms to compute) special mathematical functions, typically all related to computing probability distribution functions, often "DPQ" (Density (pdf), Probability (cdf), and Quantiles (inverse cdf)).
</p>

<p>
  Initiated and currently solely maintained by
  <a href="https://stat.ethz.ch/~maechler/">Martin Maechler</a>,
  <a href="https://www.ethz.ch/">ETH Zurich</a>,
  and member of the R Core team since primordial times.
</p>

<p>
  Traditionally, I have worked a lot in making the DPQ-functions in (base)
  R more accurate in border cases,  notably implementing
  the <tt>log=TRUE/FALSE</tt> and <tt>lower.tail = T|F</tt> possibilities.
</p>

<h3>R Packages part of the 'specfun' project on R-forge:</h3>
<!--- FIXME: rather use a nice table == make it via markdown -->
<ul>
<li>Bessel: </li>
<li>DPQ: </li>
<li>DPQmpfr: </li>
<li>dcdflib: </li>
</ul>
see the <strong>R Packages</strong> menu item on the project summary page,
currently available <a href="http://r-forge.r-project.org/R/?group_id=611">here</a>.


<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. </p>

</body>
</html>
