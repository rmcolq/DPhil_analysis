c = Channel.from(1..3)
process make_vcf {
  memory { 0.1.GB }

  input:
  val(count) from c

  output:
  file("my.vcf") into vcfs

  """
  touch my.vcf
  if [[ "${count}" == "1" ]]
  then
    echo "there" >> my.vcf
  fi
  if [[ "${count}" == "1" ]] 
  then
    echo "there" >> my.vcf
  fi
  """
}
vcfs.collectFile(name: 'all.vcf', keepHeader: true, skip: 2)
