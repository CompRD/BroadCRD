function getUrlVars()
{
    // thanks 
    // http://stackoverflow.com/questions/5422265/how-can-i-pre-populate-html-form-input-fields-from-url-parameters
    var vars = [], hash;
    var hashes = window.location.href.slice(window.location.href.indexOf('?') + 1).split('&');

    for(var i = 0; i < hashes.length; i++)
        {
         hash = hashes[i].split('=');
         vars.push(hash[0]);
         vars[hash[0]] = hash[1];
         }

     return vars;
}

function startup()
{
    args=getUrlVars();
    setSelector(urlDecode(args['sel']));
    setValue(urlDecode(args['n']));
}

function urlDecode($str) {
    if ( !$str ) return $str;
    else return decodeURIComponent(($str+'').replace(/\+/g, '%20'));
}

function keypressClear() {
    var index = document.getElementById('selvalue').value;
    if ( index != 0 ) {
	setSelector(0);
    }
}

function setSelector($index) {
    if ( !$index ) { $index=0; }
    document.getElementById('selvalue').value=$index;
    document.getElementById('selector').selectedIndex=$index;
}

function setValue($value) {
    if ( !$value ) { $value="S=11:90.26M D=1"; }
    document.getElementById('value').value=$value;
}

function setAndMaybeSubmit($new_value,$new_selvalue)
{
    var value = document.getElementById('value').value;
    var selvalue = document.getElementById('selvalue').value;

    console.log($new_value);
    console.log($new_selvalue);

    if ( $new_value != value || $new_selvalue != selvalue ) {
	setValue($new_value);
	setSelector($new_selvalue);
	if ($new_selvalue > 0 ) { submitform(); }
    }
}
function submitform()
{
    showloading()
    document.forms['subform'].submit();
}
function showloading()
{
    document.getElementById('loading').style.display='inline';
    document.getElementById('results').style.display='none';
}
function defaultform()
{
    setValue();
    setSelector();
    submitform();
}
function scalehack()
{
    var height = 0.5*document.getElementById('imageout').height;
    var width = 0.5*document.getElementById('imageout').width;
    var adj = 1.0;
    if ( height < 90 ) { adj = 90.0 / height; } else if ( width < 90 ) { adj = 90.0 / width; } else { adj = 1.0; }
    height *= adj;
    width *= adj;
    document.getElementById('imageout').height = height;
    document.getElementById('imageout').width = width;
}
