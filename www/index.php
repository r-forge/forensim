<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">
<head>
<!-- This is the project specific website template --><!-- It can be changed as liked or replaced by other content --><?php $domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>
  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
  <title>&lt;?php echo $group_name; ?&gt;</title>
  <link href="%3C?php%20echo%20$themeroot;%20?%3Estyles/estilo1.css"
 rel="stylesheet" type="text/css" />
</head>
<body>
';
?&gt;<!-- R-Forge Logo -->
<table border="0" cellpadding="0" cellspacing="0" width="100%">
  <tbody>
    <tr>
      <td><a href="/"><img
 src="%3C?php%20echo%20$themeroot;%20?%3E/images/logo.png"
 alt="R-Forge Logo" border="0" /> </a> </td>
    </tr>
  </tbody>
</table>
<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like --><?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>
<!-- end of project description -->
<p> No content added.&nbsp; </p>
<p> The <strong>project summary page</strong> you can find <a
 href="http://%3C?php%20echo%20$domain;%20?%3E/projects/%3C?php%20echo%20$group_name;%20?%3E/"><strong>here</strong></a>.
</p>
</body>
</html>
